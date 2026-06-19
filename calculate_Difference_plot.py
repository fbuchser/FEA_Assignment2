import re
import numpy as np
import pyvista as pv
from pathlib import Path
from scipy.io import loadmat


# ── Paths ─────────────────────────────────────────────────────────────────────
INPUT_PATH = Path(r"data\FemurInput.mat")
OWN_RESULTS_PATH = Path(r"results\sideways_fall_results.mat")
ANSYS_EXPORT_PATH = Path(r"results_ansys\numerical")


ANSYS_FILES = {
    "disp": ANSYS_EXPORT_PATH / "displacements.lis",
    "stress": ANSYS_EXPORT_PATH / "prin_stress.lis",
    "strain": ANSYS_EXPORT_PATH / "prin_strain.lis",
}


# ── Plot settings ─────────────────────────────────────────────────────────────
PLOT_CMAP = "coolwarm"
N_LABELS = 5
SHOW_EDGES = False
TITLE_STR = "Sideways Fall"


# ── Regex ─────────────────────────────────────────────────────────────────────
NUM_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")
ELEMENT_RE = re.compile(r"ELEMENT=\s*(\d+)")
NODE_ID_RE = re.compile(r"^\s*(\d+)\s+(.+)$")


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
    elem_ids = []
    values = []
    current_eid = None
    current_rows = []

    def flush():
        if current_eid is not None and current_rows:
            arr = np.array(current_rows, dtype=float)
            elem_ids.append(current_eid)
            values.append(arr.mean(axis=0))

    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m_elem = ELEMENT_RE.search(line)
            if m_elem:
                flush()
                current_eid = int(m_elem.group(1))
                current_rows = []
                continue

            stripped = line.strip()
            if not stripped or stripped.startswith("NODE") or stripped.startswith("*") \
               or stripped.startswith("PRINT") or stripped.startswith("LOAD") \
               or stripped.startswith("TIME") or stripped.startswith("THE") \
               or stripped.startswith("POST"):
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


def align_by_ids_strict(target_ids, source_ids, source_values, field_name):
    """
    Strict alignment: every target ID must exist in source.
    Use for element-wise fields (stress, strain) where coverage must be complete.
    """
    src_map = {int(nid): i for i, nid in enumerate(source_ids)}
    missing = [nid for nid in target_ids if int(nid) not in src_map]
    if missing:
        raise ValueError(
            f"Missing {len(missing)} IDs in '{field_name}'. First few: {missing[:10]}"
        )
    idx = np.array([src_map[int(nid)] for nid in target_ids], dtype=int)
    return source_values[idx]


def align_by_ids_partial(target_ids, source_ids, source_values,
                         field_name="field", fill_value=0.0):
    """
    Partial alignment: missing source IDs are filled with fill_value (default 0.0).
    Use for nodal displacements when ANSYS only exports a subset of nodes
    (e.g. only non-zero / non-constrained nodes).

    Returns:
        aligned_values : np.ndarray, same length as target_ids
        missing_ids    : list of IDs in target_ids not found in source_ids
    """
    target_ids = np.asarray(target_ids).astype(int).ravel()
    source_ids = np.asarray(source_ids).astype(int).ravel()
    source_values = np.asarray(source_values)

    out_shape = (len(target_ids),) + source_values.shape[1:]
    aligned = np.full(out_shape, fill_value, dtype=float)

    src_map = {int(nid): i for i, nid in enumerate(source_ids)}
    missing_ids = []
    matched = 0

    for i, nid in enumerate(target_ids):
        j = src_map.get(int(nid))
        if j is not None:
            aligned[i] = source_values[j]
            matched += 1
        else:
            missing_ids.append(int(nid))

    missing_pct = 100.0 * len(missing_ids) / max(len(target_ids), 1)
    print(
        f"  {field_name}: matched {matched}/{len(target_ids)} IDs, "
        f"missing {len(missing_ids)} ({missing_pct:.1f}%) → filled with {fill_value}"
    )

    if missing_pct > 20.0:
        print(
            f"  WARNING: More than 20% of {field_name} nodes are missing from ANSYS output.\n"
            f"  Difference plots for displacement may be unreliable.\n"
            f"  Consider re-exporting ANSYS displacements for ALL nodes (NSEL,ALL before PRNSOL,U,COMP)."
        )

    return aligned, missing_ids


def align_with_nan(target_ids, source_ids, source_values, field_name):
    src_map = {int(eid): i for i, eid in enumerate(source_ids)}
    n_cols = source_values.shape[1]
    out = np.full((len(target_ids), n_cols), np.nan, dtype=float)
    found_mask = np.array([int(eid) in src_map for eid in target_ids])
    for i, eid in enumerate(target_ids):
        if found_mask[i]:
            out[i] = source_values[src_map[int(eid)]]
    missing_idx = np.where(~found_mask)[0]
    if len(missing_idx):
        print(f"  Warning: {field_name}: matched {found_mask.sum()}/{len(target_ids)}; {len(missing_idx)} left as NaN.")
    return out


def build_mesh(coords, conn):
    n_elem = conn.shape[0]
    cells = np.hstack([np.full((n_elem, 1), 4, dtype=np.int64), conn]).ravel()
    cell_types = np.full(n_elem, pv.CellType.TETRA, dtype=np.uint8)
    return pv.UnstructuredGrid(cells, cell_types, coords)


def node_to_element_average(node_vals, conn):
    return np.asarray(node_vals)[conn].mean(axis=1)


def symmetric_limits(*arrays):
    valid = []
    for a in arrays:
        if a is None:
            continue
        arr = np.asarray(a).ravel()
        arr = arr[np.isfinite(arr)]
        if arr.size:
            valid.append(arr)
    if not valid:
        return -1.0, 1.0
    vmax = max(np.nanmax(np.abs(a)) for a in valid)
    return (-vmax, vmax) if vmax != 0 else (-1.0, 1.0)


def add_diff_plot(plotter, mesh, values, title, scalar_name, clim):
    m = mesh.copy()
    m.cell_data[scalar_name] = values
    plotter.add_mesh(
        m,
        scalars=scalar_name,
        cmap=PLOT_CMAP,
        clim=clim,
        show_edges=SHOW_EDGES,
        scalar_bar_args={"title": scalar_name, "color": "black", "n_labels": N_LABELS},
    )
    plotter.add_text(title, font_size=10, color="black")
    plotter.add_axes(xlabel="X [mm]", ylabel="Y [mm]", zlabel="Z [mm]")


def plot_difference_fields(mesh, diff_fields, title_str):
    disp_clim = symmetric_limits(diff_fields["du_mag"])
    stress_clim = symmetric_limits(diff_fields["ds1"], diff_fields["ds3"])
    strain_clim = symmetric_limits(diff_fields["de1"], diff_fields["de3"])
    uz_clim = symmetric_limits(diff_fields["duz"])

    pl = pv.Plotter(shape=(2, 3), window_size=(3000, 1600))
    pl.set_background("white")

    panels = [
        (0, 0, diff_fields["du_mag"], "Displacement magnitude relative error", "Δ||u|| [%]", disp_clim),
        (0, 1, diff_fields["ds1"], "Max principal stress relative error", "Δs1 [%]", stress_clim),
        (0, 2, diff_fields["ds3"], "Min principal stress relative error", "Δs3 [%]", stress_clim),
        (1, 0, diff_fields["de1"], "Max principal strain relative error", "Δe1 [%]", strain_clim),
        (1, 1, diff_fields["de3"], "Min principal strain relative error", "Δe3 [%]", strain_clim),
        (1, 2, diff_fields["duz"], "Z-displacement relative error", "Δuz [%]", uz_clim),
    ]

    for row, col, values, title, scalar_name, clim in panels:
        pl.subplot(row, col)
        add_diff_plot(pl, mesh, values, title, scalar_name, clim)

    pl.add_text(
        f'ANSYS vs own solver relative errors [%] – "{title_str}"\n(defined as (own - ANSYS) / |ANSYS| × 100)',
        font_size=13,
        color="black",
        position="upper_edge",
    )
    pl.show()


def compute_scalar_metrics(own_vals, ref_vals, name, unit=""):
    own_vals = np.asarray(own_vals, dtype=float).ravel()
    ref_vals = np.asarray(ref_vals, dtype=float).ravel()
    valid = np.isfinite(own_vals) & np.isfinite(ref_vals)

    if not np.any(valid):
        return {
            "name": name,
            "n": 0,
            "max_abs": np.nan,
            "l2_abs": np.nan,
            "rel_l2": np.nan,
            "rmse": np.nan,
            "rel_rmse": np.nan,
            "unit": unit,
        }

    diff = own_vals[valid] - ref_vals[valid]
    ref = ref_vals[valid]

    max_abs = np.max(np.abs(diff))
    l2_abs = np.linalg.norm(diff)
    ref_l2 = np.linalg.norm(ref)
    rmse = np.sqrt(np.mean(diff**2))
    ref_rms = np.sqrt(np.mean(ref**2))

    return {
        "name": name,
        "n": int(np.sum(valid)),
        "max_abs": max_abs,
        "ref_max": np.max(np.abs(ref)),
        "l2_abs": l2_abs,
        "rel_l2": l2_abs / ref_l2 if ref_l2 != 0 else np.nan,
        "rmse": rmse,
        "rel_rmse": rmse / ref_rms if ref_rms != 0 else np.nan,
        "unit": unit,
    }


def compute_vector_metrics(own_vec, ref_vec, name, unit=""):
    own_vec = np.asarray(own_vec, dtype=float)
    ref_vec = np.asarray(ref_vec, dtype=float)
    valid = np.all(np.isfinite(own_vec), axis=1) & np.all(np.isfinite(ref_vec), axis=1)

    if not np.any(valid):
        return {
            "name": name,
            "n": 0,
            "max_abs": np.nan,
            "l2_abs": np.nan,
            "rel_l2": np.nan,
            "rmse": np.nan,
            "rel_rmse": np.nan,
            "unit": unit,
        }

    diff = own_vec[valid] - ref_vec[valid]
    ref = ref_vec[valid]

    diff_flat = diff.ravel()
    ref_flat = ref.ravel()
    diff_mag = np.linalg.norm(diff, axis=1)

    max_abs = np.max(diff_mag)
    l2_abs = np.linalg.norm(diff_flat)
    ref_l2 = np.linalg.norm(ref_flat)
    rmse = np.sqrt(np.mean(diff_flat**2))
    ref_rms = np.sqrt(np.mean(ref_flat**2))

    return {
        "name": name,
        "n": int(np.sum(valid)),
        "max_abs": max_abs,
        "ref_max": np.max(np.abs(ref)),
        "l2_abs": l2_abs,
        "rel_l2": l2_abs / ref_l2 if ref_l2 != 0 else np.nan,
        "rmse": rmse,
        "rel_rmse": rmse / ref_rms if ref_rms != 0 else np.nan,
        "unit": unit,
    }


def print_metrics_table(metrics_list):
    print("Numerical error metrics (own solver vs ANSYS)")
    print("-" * 110)
    header = f"{'Quantity':<24} {'n':>8} {'Max abs err':>16} {'L2 abs err':>16} {'Rel L2 [%]':>12} {'RMSE':>16} {'Rel RMSE [%]':>14}"
    print(header)
    print("-" * 110)

    for m in metrics_list:
        unit = f" {m['unit']}" if m['unit'] else ""
        max_abs = f"{m['max_abs']:.6g}{unit}" if np.isfinite(m['max_abs']) else "nan"
        l2_abs = f"{m['l2_abs']:.6g}{unit}" if np.isfinite(m['l2_abs']) else "nan"
        rmse = f"{m['rmse']:.6g}{unit}" if np.isfinite(m['rmse']) else "nan"
        rel_l2 = f"{100*m['rel_l2']:.4f}" if np.isfinite(m['rel_l2']) else "nan"
        rel_rmse = f"{100*m['rel_rmse']:.4f}" if np.isfinite(m['rel_rmse']) else "nan"

        print(f"{m['name']:<24} {m['n']:>8d} {max_abs:>16} {l2_abs:>16} {rel_l2:>12} {rmse:>16} {rel_rmse:>14}")

    print("-" * 110)


# ── Main ───────────────────────────────────────────────────────────────────────
inp = loadmat(INPUT_PATH)
nodes = inp["nodes"]
elements = inp["elements"]

coords_solver = nodes[:, 1:4].astype(float)
conn = (elements[:, 1:5] - 1).astype(int)
node_ids_solver = nodes[:, 0].astype(int)
elem_ids_solver = elements[:, 0].astype(int)
n_nodes = coords_solver.shape[0]
n_elem = conn.shape[0]
mesh = build_mesh(coords_solver, conn)

print(f"Mesh: {n_nodes} nodes, {n_elem} elements")

own = loadmat(OWN_RESULTS_PATH)
U_own = read_first_available(own, "U", "displacement", "disp").reshape(n_nodes, 3)
s1_own = read_first_available(own, "sigma1", "s1")
s3_own = read_first_available(own, "sigma3", "s3")
e1_own = read_first_available(own, "epsi1", "epsilon1", "strain1", "e1")
e3_own = read_first_available(own, "epsi3", "epsilon3", "strain3", "e3")

print("Loading ANSYS files...")
disp_ids_ansys, U_ansys = load_apdl_displacements(ANSYS_FILES["disp"])
stress_ids_ansys, stress_prin = load_element_nodal_listing(ANSYS_FILES["stress"])
strain_ids_ansys, strain_prin = load_element_nodal_listing(ANSYS_FILES["strain"])

print(f"ANSYS disp nodes: {len(disp_ids_ansys)}, stress elems: {len(stress_ids_ansys)}, strain elems: {len(strain_ids_ansys)}")

# Displacements: ANSYS may only export a subset of nodes (e.g. non-constrained).
# Missing nodes are filled with 0.0 (assumed zero displacement at BCs).
U_ansys_aligned, _missing_disp = align_by_ids_partial(
    node_ids_solver, disp_ids_ansys, U_ansys,
    field_name="displacements", fill_value=0.0
)

# Stress / strain: use NaN-tolerant alignment (same behaviour as before).
stress_ansys_aligned = align_with_nan(elem_ids_solver, stress_ids_ansys, stress_prin, "principal stresses")
strain_ansys_aligned = align_with_nan(elem_ids_solver, strain_ids_ansys, strain_prin, "principal strains")

u_mag_own = np.linalg.norm(U_own, axis=1)
u_mag_ansys = np.linalg.norm(U_ansys_aligned, axis=1)

s1_ansys = stress_ansys_aligned[:, 0]
s3_ansys = stress_ansys_aligned[:, 2]
e1_ansys = strain_ansys_aligned[:, 0]
e3_ansys = strain_ansys_aligned[:, 2]


def rel_error_pct(own, ref, eps_frac=1e-12):
    """Relative error in %, |own - ref| / (|ref| + eps) * 100.
    eps guards against division by zero at near-zero reference values."""
    own = np.asarray(own, dtype=float)
    ref = np.asarray(ref, dtype=float)
    eps = eps_frac * (np.nanmax(np.abs(ref)) if np.any(np.isfinite(ref)) else 1.0)
    return (own - ref) / (np.abs(ref) + eps) * 100.0


# Nodal relative errors → averaged onto elements for plotting
du_mag_elem  = node_to_element_average(rel_error_pct(u_mag_own, u_mag_ansys), conn)
duz_elem     = node_to_element_average(rel_error_pct(U_own[:, 2], U_ansys_aligned[:, 2]), conn)

ds1_elem = rel_error_pct(s1_own, s1_ansys)
ds3_elem = rel_error_pct(s3_own, s3_ansys)
de1_elem = rel_error_pct(e1_own, e1_ansys)
de3_elem = rel_error_pct(e3_own, e3_ansys)

diff_fields = {
    "du_mag": du_mag_elem,
    "duz": duz_elem,
    "ds1": ds1_elem,
    "ds3": ds3_elem,
    "de1": de1_elem,
    "de3": de3_elem,
}

metrics = [
    compute_vector_metrics(U_own, U_ansys_aligned, "Displacement vector", unit="mm"),
    compute_scalar_metrics(np.linalg.norm(U_own, axis=1), np.linalg.norm(U_ansys_aligned, axis=1), "Displacement magnitude", unit="mm"),
    compute_scalar_metrics(U_own[:, 2], U_ansys_aligned[:, 2], "Z displacement", unit="mm"),
    compute_scalar_metrics(s1_own, s1_ansys, "Max principal stress", unit="MPa"),
    compute_scalar_metrics(s3_own, s3_ansys, "Min principal stress", unit="MPa"),
    compute_scalar_metrics(e1_own, e1_ansys, "Max principal strain"),
    compute_scalar_metrics(e3_own, e3_ansys, "Min principal strain"),
]

print_metrics_table(metrics)
plot_difference_fields(mesh, diff_fields, TITLE_STR)

# Check if absolute stress errors are actually large, or just inflated by near-zero reference
valid = np.isfinite(s1_ansys) & np.isfinite(s1_own)

abs_err_s1 = np.abs(s1_own[valid] - s1_ansys[valid])
ref_s1 = np.abs(s1_ansys[valid])

# What fraction of elements have BOTH large relative AND large absolute error?
large_rel  = (abs_err_s1 / (ref_s1 + 1e-12 * ref_s1.max())) > 0.20  # >20% rel error
large_abs  = abs_err_s1 > 1.0  # >1 MPa absolute error

print("Elements with >20% relative error:", large_rel.sum(), "/", valid.sum())
print("Elements with >1 MPa absolute error:", large_abs.sum(), "/", valid.sum())
print("Both large rel AND large abs:", (large_rel & large_abs).sum())
print()
print("Max absolute stress error [MPa]:", abs_err_s1.max())
print("Median absolute stress error [MPa]:", np.median(abs_err_s1))
print("90th percentile absolute stress error [MPa]:", np.percentile(abs_err_s1, 90))
print()
# Check if near-zero reference values are inflating the percentage
print("Median |ANSYS s1| [MPa]:", np.median(ref_s1))
print("Fraction of elements where |ANSYS s1| < 1 MPa:", (ref_s1 < 1.0).mean())