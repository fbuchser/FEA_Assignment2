import re
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors
from pathlib import Path
from scipy.io import loadmat


# ── LaTeX / font settings ─────────────────────────────────────────────────────
mpl.rcParams.update({
    "text.usetex":         False,
    "font.family":         "serif",
    "mathtext.fontset":    "cm",
    "axes.labelsize":      10,
    "xtick.labelsize":     8,
    "ytick.labelsize":     8,
    "figure.facecolor":    "white",
    "axes.facecolor":      "white",
})


# ── Paths ─────────────────────────────────────────────────────────────────────
INPUT_PATH       = Path(r"data\FemurInput.mat")
OWN_RESULTS_PATH = Path(r"results\sideways_fall_results.mat")
ANSYS_EXPORT_PATH= Path(r"results_ansys\numerical")

ANSYS_FILES = {
    "disp":   ANSYS_EXPORT_PATH / "displacements.lis",
    "stress": ANSYS_EXPORT_PATH / "prin_stress.lis",
    "strain": ANSYS_EXPORT_PATH / "prin_strain.lis",
}


# ── Plot settings ─────────────────────────────────────────────────────────────
PLOT_CMAP  = "coolwarm"
SHOW_EDGES = True
TITLE_STR  = "Sideways Fall"
OUTPUT_DIR = Path(r"results_ansys\compare")


# ── Regex ─────────────────────────────────────────────────────────────────────
NUM_RE     = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")
ELEMENT_RE = re.compile(r"ELEMENT=\s*(\d+)")
NODE_ID_RE = re.compile(r"^\s*(\d+)\s+(.+)$")


MATPLOTLIB_VIEWS = [
    (15, 145, "Anterior"),
    (15,  -5, "Lateral-Posterior"),
]


# Panel definitions: (field_key, colorbar_label_latex, title_latex, safe_filename)
PANELS = [
    ("du_mag", r"$\Delta\Vert\mathbf{u}\Vert\ [\%]$",
               r"Displacement magnitude relative error --- " + TITLE_STR,
               "displacement_magnitude_relative_error"),
    ("ds1",    r"$\Delta\sigma_1\ \mathrm{[MPa]}$",
               r"Max. principal stress absolute difference --- " + TITLE_STR,
               "max_principal_stress_relative_error"),
    ("ds3",    r"$\Delta\sigma_3\ \mathrm{[MPa]}$",
               r"Min. principal stress absolute difference --- " + TITLE_STR,
               "min_principal_stress_relative_error"),
    ("de1",    r"$\Delta\varepsilon_1\ [-]$",
               r"Max. principal strain absolute difference --- " + TITLE_STR,
               "max_principal_strain_relative_error"),
    ("de3",    r"$\Delta\varepsilon_3\ [-]$",
               r"Min. principal strain absolute difference --- " + TITLE_STR,
               "min_principal_strain_relative_error"),
    ("duz",    r"$\Delta u_z\ \mathrm{[mm]}$",
               r"Z-displacement absolute difference --- " + TITLE_STR,
               "z_displacement_relative_error"),
]


# ── Helpers ───────────────────────────────────────────────────────────────────
def _to_float(tok):
    return float(tok.replace("D", "E").replace("d", "E"))

def parse_any_floats(s):
    return [_to_float(t) for t in NUM_RE.findall(s)]

def read_first_available(matdict, *keys):
    for key in keys:
        if key in matdict:
            return np.asarray(matdict[key]).ravel()
    raise KeyError(f"None of these keys were found: {keys}")

def extract_listing_rows(file_path, n_value_cols):
    rows = []
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            vals = parse_any_floats(s)
            if len(vals) < 1 + n_value_cols:
                continue
            first = vals[0]
            if abs(first - round(first)) > 1e-9:
                continue
            row = [int(round(first))] + vals[1:1 + n_value_cols]
            rows.append(row)
    if not rows:
        raise ValueError(f"No numeric rows found in {file_path}.")
    arr = np.asarray(rows, dtype=float)
    ids = arr[:, 0].astype(int)
    _, unique_idx = np.unique(ids, return_index=True)
    return arr[np.sort(unique_idx)]

def load_apdl_displacements(file_path):
    data = extract_listing_rows(file_path, n_value_cols=3)
    return data[:, 0].astype(int), data[:, 1:4]

def load_element_nodal_listing(file_path):
    elem_ids, values = [], []
    current_eid, current_rows = None, []

    def flush():
        if current_eid is not None and current_rows:
            elem_ids.append(current_eid)
            values.append(np.array(current_rows, dtype=float).mean(axis=0))

    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m_elem = ELEMENT_RE.search(line)
            if m_elem:
                flush()
                current_eid, current_rows = int(m_elem.group(1)), []
                continue
            stripped = line.strip()
            if not stripped or stripped.startswith(("NODE","*","PRINT","LOAD","TIME","THE","POST")):
                continue
            m_node = NODE_ID_RE.match(line)
            if m_node and current_eid is not None:
                nums = parse_any_floats(m_node.group(2))
                if len(nums) >= 3:
                    current_rows.append(nums[:3])
    flush()
    if not elem_ids:
        raise ValueError(f"No element blocks found in {file_path}")
    print(f"  Parsed {len(elem_ids)} elements from {file_path.name}")
    return np.array(elem_ids, dtype=int), np.array(values, dtype=float)

def align_by_ids_partial(target_ids, source_ids, source_values,
                         field_name="field", fill_value=0.0):
    target_ids    = np.asarray(target_ids).astype(int).ravel()
    source_ids    = np.asarray(source_ids).astype(int).ravel()
    source_values = np.asarray(source_values)
    out_shape     = (len(target_ids),) + source_values.shape[1:]
    aligned       = np.full(out_shape, fill_value, dtype=float)
    src_map       = {int(nid): i for i, nid in enumerate(source_ids)}
    missing_ids, matched = [], 0
    for i, nid in enumerate(target_ids):
        j = src_map.get(int(nid))
        if j is not None:
            aligned[i] = source_values[j]; matched += 1
        else:
            missing_ids.append(int(nid))
    missing_pct = 100.0 * len(missing_ids) / max(len(target_ids), 1)
    print(f"  {field_name}: matched {matched}/{len(target_ids)}, "
          f"missing {len(missing_ids)} ({missing_pct:.1f}%) -> filled with {fill_value}")
    if missing_pct > 20.0:
        print(f"  WARNING: >20% of {field_name} nodes missing from ANSYS output.")
    return aligned, missing_ids

def align_with_nan(target_ids, source_ids, source_values, field_name):
    src_map    = {int(eid): i for i, eid in enumerate(source_ids)}
    n_cols     = source_values.shape[1]
    out        = np.full((len(target_ids), n_cols), np.nan, dtype=float)
    found_mask = np.array([int(eid) in src_map for eid in target_ids])
    for i, eid in enumerate(target_ids):
        if found_mask[i]:
            out[i] = source_values[src_map[int(eid)]]
    missing = np.sum(~found_mask)
    if missing:
        print(f"  Warning: {field_name}: matched {found_mask.sum()}/{len(target_ids)}; {missing} left as NaN.")
    return out

def build_mesh(coords, conn):
    n_elem     = conn.shape[0]
    cells      = np.hstack([np.full((n_elem, 1), 4, dtype=np.int64), conn]).ravel()
    cell_types = np.full(n_elem, pv.CellType.TETRA, dtype=np.uint8)
    return pv.UnstructuredGrid(cells, cell_types, coords)

def node_to_element_average(node_vals, conn):
    return np.asarray(node_vals)[conn].mean(axis=1)

def rel_error_pct(own, ref, eps_frac=1e-12):
    own = np.asarray(own, dtype=float)
    ref = np.asarray(ref, dtype=float)
    eps = eps_frac * (np.nanmax(np.abs(ref)) if np.any(np.isfinite(ref)) else 1.0)
    return (own - ref) / (np.abs(ref) + eps) * 100.0

def abs_diff(own, ref):
    """Absolute difference: own solver minus ANSYS."""
    return np.asarray(own, dtype=float) - np.asarray(ref, dtype=float)


# ── Plotting ──────────────────────────────────────────────────────────────────
def plot_difference_fields(mesh, diff_fields, title_str):
    OUTPUT_DIR.mkdir(exist_ok=True)

    points = np.asarray(mesh.points)
    cells  = mesh.cells.reshape(-1, 5)[:, 1:5]

    tri_faces  = np.vstack([cells[:, [0,1,2]], cells[:, [0,1,3]],
                             cells[:, [0,2,3]], cells[:, [1,2,3]]])
    tri_sorted = np.sort(tri_faces, axis=1)
    uniq, counts = np.unique(tri_sorted, axis=0, return_counts=True)
    boundary_set  = {tuple(r) for r in uniq[counts == 1]}
    boundary_faces= np.array(
        [r for r in tri_faces if tuple(sorted(r)) in boundary_set], dtype=int)

    x, y, z = points[:, 0], points[:, 1], points[:, 2]

    for field_key, cbar_label, suptitle, safe_name in PANELS:
        elem_values = np.asarray(diff_fields[field_key], dtype=float)

        node_acc = np.zeros(len(points)); node_cnt = np.zeros(len(points))
        for e, tet in enumerate(cells):
            node_acc[tet] += elem_values[e]; node_cnt[tet] += 1.0
        node_values = node_acc / np.maximum(node_cnt, 1.0)

        finite = np.isfinite(node_values)
        vmax   = np.nanmax(np.abs(node_values[finite])) if np.any(finite) else 1.0
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0
        norm = colors.Normalize(vmin=-vmax, vmax=vmax)
        cmap = plt.get_cmap(PLOT_CMAP)

        fig  = plt.figure(figsize=(10.5, 5.4), facecolor="white")
        axes = [fig.add_subplot(1, 2, i + 1, projection="3d") for i in range(2)]

        tri_vals   = node_values[boundary_faces].mean(axis=1)
        tri_colors = cmap(norm(tri_vals))

        for ax, (elev, azim, view_label) in zip(axes, MATPLOTLIB_VIEWS):
            surf = ax.plot_trisurf(
                x, y, z,
                triangles=boundary_faces,
                linewidth=0.15 if SHOW_EDGES else 0.0,
                edgecolor="grey" if SHOW_EDGES else "none",
                antialiased=True,
                shade=False,
            )
            surf.set_facecolors(tri_colors)

            ax.view_init(elev=elev, azim=azim)
            ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))

            for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
                pane.fill = False
                pane.set_edgecolor((0.85, 0.85, 0.85, 1.0))
            ax.set_facecolor("white")
            ax.grid(False)

            # Axis labels — tight labelpad
            ax.set_xlabel(r"$X\ \mathrm{[mm]}$", labelpad=1)
            ax.set_ylabel(r"$Y\ \mathrm{[mm]}$", labelpad=1)
            ax.set_zlabel(r"$Z\ \mathrm{[mm]}$", labelpad=1)

            # Tick spacing
            for arr, setter in [(x, ax.set_xticks), (y, ax.set_yticks), (z, ax.set_zticks)]:
                t0 = 20.0 * np.floor(np.nanmin(arr) / 20.0)
                t1 = 20.0 * np.ceil(np.nanmax(arr)  / 20.0)
                setter(np.arange(t0, t1 + 1e-9, 20.0))

            # Hide last X tick label on left subplot (removes overhanging "100")
            if azim == 145:
                ax.xaxis.get_major_ticks()[-1].label1.set_visible(False)

            # Hide last X tick label on right subplot too (the stray "100")
            if azim == -5:
                ax.xaxis.get_major_ticks()[-1].label1.set_visible(False)

            # Reduced pad pulls tick labels closer to the axis
            ax.tick_params(axis="x", labelsize=8, pad=-2)
            ax.tick_params(axis="y", labelsize=8, pad=-2)
            ax.tick_params(axis="z", labelsize=8, pad=-2)

        fig.subplots_adjust(left=0.06, right=0.78, bottom=0.08, top=0.88, wspace=-0.20)

        cax      = fig.add_axes([0.735, 0.18, 0.014, 0.58])
        mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        mappable.set_array(node_values)
        cbar     = fig.colorbar(mappable, cax=cax)
        cbar.set_label(cbar_label, rotation=90, labelpad=10, fontsize=10)
        cbar.ax.tick_params(labelsize=8)

        fig.suptitle(suptitle, fontsize=11, y=0.96)

        out_path = OUTPUT_DIR / f"{safe_name}.png"
        fig.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.3, facecolor="white")
        plt.close(fig)
        print(f"Saved: {out_path}")


# ── Metrics ───────────────────────────────────────────────────────────────────
def compute_scalar_metrics(own_vals, ref_vals, name, unit=""):
    own_vals = np.asarray(own_vals, dtype=float).ravel()
    ref_vals = np.asarray(ref_vals, dtype=float).ravel()
    valid    = np.isfinite(own_vals) & np.isfinite(ref_vals)
    if not np.any(valid):
        return {"name": name, "n": 0, "max_abs": np.nan, "l2_abs": np.nan,
                "rel_l2": np.nan, "rmse": np.nan, "rel_rmse": np.nan, "unit": unit}
    diff = own_vals[valid] - ref_vals[valid]; ref = ref_vals[valid]
    l2d  = np.linalg.norm(diff); l2r = np.linalg.norm(ref)
    rms  = np.sqrt(np.mean(diff**2)); rmsr= np.sqrt(np.mean(ref**2))
    return {"name": name, "n": int(np.sum(valid)),
            "max_abs": np.max(np.abs(diff)), "ref_max": np.max(np.abs(ref)),
            "l2_abs": l2d, "rel_l2": l2d/l2r if l2r else np.nan,
            "rmse": rms,   "rel_rmse": rms/rmsr if rmsr else np.nan, "unit": unit}

def compute_vector_metrics(own_vec, ref_vec, name, unit=""):
    own_vec = np.asarray(own_vec, dtype=float)
    ref_vec = np.asarray(ref_vec, dtype=float)
    valid   = np.all(np.isfinite(own_vec), axis=1) & np.all(np.isfinite(ref_vec), axis=1)
    if not np.any(valid):
        return {"name": name, "n": 0, "max_abs": np.nan, "l2_abs": np.nan,
                "rel_l2": np.nan, "rmse": np.nan, "rel_rmse": np.nan, "unit": unit}
    diff = own_vec[valid] - ref_vec[valid]; ref = ref_vec[valid]
    df   = diff.ravel(); rf = ref.ravel()
    l2d  = np.linalg.norm(df); l2r = np.linalg.norm(rf)
    rms  = np.sqrt(np.mean(df**2)); rmsr= np.sqrt(np.mean(rf**2))
    return {"name": name, "n": int(np.sum(valid)),
            "max_abs": np.max(np.linalg.norm(diff, axis=1)), "ref_max": np.max(np.abs(ref)),
            "l2_abs": l2d, "rel_l2": l2d/l2r if l2r else np.nan,
            "rmse": rms,   "rel_rmse": rms/rmsr if rmsr else np.nan, "unit": unit}

def print_metrics_table(metrics_list):
    print("\nNumerical error metrics (own solver vs ANSYS)")
    print("-" * 110)
    print(f"{'Quantity':<24} {'n':>8} {'Max abs err':>16} {'L2 abs err':>16} {'Rel L2 [%]':>12} {'RMSE':>16} {'Rel RMSE [%]':>14}")
    print("-" * 110)
    for m in metrics_list:
        unit     = f" {m['unit']}" if m['unit'] else ""
        max_abs  = f"{m['max_abs']:.6g}{unit}" if np.isfinite(m['max_abs']) else "nan"
        l2_abs   = f"{m['l2_abs']:.6g}{unit}"  if np.isfinite(m['l2_abs'])  else "nan"
        rmse     = f"{m['rmse']:.6g}{unit}"     if np.isfinite(m['rmse'])    else "nan"
        rel_l2   = f"{100*m['rel_l2']:.4f}"    if np.isfinite(m['rel_l2'])  else "nan"
        rel_rmse = f"{100*m['rel_rmse']:.4f}"  if np.isfinite(m['rel_rmse'])else "nan"
        print(f"{m['name']:<24} {m['n']:>8d} {max_abs:>16} {l2_abs:>16} {rel_l2:>12} {rmse:>16} {rel_rmse:>14}")
    print("-" * 110)


# ── Main ───────────────────────────────────────────────────────────────────────
inp      = loadmat(INPUT_PATH)
nodes    = inp["nodes"]
elements = inp["elements"]

coords_solver   = nodes[:, 1:4].astype(float)
conn            = (elements[:, 1:5] - 1).astype(int)
node_ids_solver = nodes[:, 0].astype(int)
elem_ids_solver = elements[:, 0].astype(int)
n_nodes = coords_solver.shape[0]; n_elem = conn.shape[0]
mesh    = build_mesh(coords_solver, conn)
print(f"Mesh: {n_nodes} nodes, {n_elem} elements")

own   = loadmat(OWN_RESULTS_PATH)
U_own = read_first_available(own, "U", "displacement", "disp").reshape(n_nodes, 3)
s1_own= read_first_available(own, "sigma1", "s1")
s3_own= read_first_available(own, "sigma3", "s3")
e1_own= read_first_available(own, "epsi1", "epsilon1", "strain1", "e1")
e3_own= read_first_available(own, "epsi3", "epsilon3", "strain3", "e3")

print("Loading ANSYS files...")
disp_ids_ansys,   U_ansys       = load_apdl_displacements(ANSYS_FILES["disp"])
stress_ids_ansys, stress_prin   = load_element_nodal_listing(ANSYS_FILES["stress"])
strain_ids_ansys, strain_prin   = load_element_nodal_listing(ANSYS_FILES["strain"])
print(f"ANSYS disp nodes: {len(disp_ids_ansys)}, stress elems: {len(stress_ids_ansys)}, strain elems: {len(strain_ids_ansys)}")

U_ansys_aligned, _ = align_by_ids_partial(
    node_ids_solver, disp_ids_ansys, U_ansys, field_name="displacements", fill_value=0.0)
stress_ansys_aligned = align_with_nan(elem_ids_solver, stress_ids_ansys, stress_prin, "principal stresses")
strain_ansys_aligned = align_with_nan(elem_ids_solver, strain_ids_ansys, strain_prin, "principal strains")

u_mag_own   = np.linalg.norm(U_own,           axis=1)
u_mag_ansys = np.linalg.norm(U_ansys_aligned, axis=1)
s1_ansys = stress_ansys_aligned[:, 0]; s3_ansys = stress_ansys_aligned[:, 2]
e1_ansys = strain_ansys_aligned[:, 0]; e3_ansys = strain_ansys_aligned[:, 2]

diff_fields = {
    "du_mag": node_to_element_average(rel_error_pct(u_mag_own, u_mag_ansys), conn),
    "duz":    node_to_element_average(abs_diff(U_own[:, 2], U_ansys_aligned[:, 2]), conn),
    "ds1":    abs_diff(s1_own, s1_ansys),
    "ds3":    abs_diff(s3_own, s3_ansys),
    "de1":    abs_diff(e1_own, e1_ansys),
    "de3":    abs_diff(e3_own, e3_ansys),
}

metrics = [
    compute_vector_metrics(U_own, U_ansys_aligned,           "Displacement vector",   unit="mm"),
    compute_scalar_metrics(u_mag_own, u_mag_ansys,           "Displacement magnitude",unit="mm"),
    compute_scalar_metrics(U_own[:,2], U_ansys_aligned[:,2], "Z displacement",        unit="mm"),
    compute_scalar_metrics(s1_own, s1_ansys,                 "Max principal stress",  unit="MPa"),
    compute_scalar_metrics(s3_own, s3_ansys,                 "Min principal stress",  unit="MPa"),
    compute_scalar_metrics(e1_own, e1_ansys,                 "Max principal strain"),
    compute_scalar_metrics(e3_own, e3_ansys,                 "Min principal strain"),
]

print_metrics_table(metrics)
plot_difference_fields(mesh, diff_fields, TITLE_STR)