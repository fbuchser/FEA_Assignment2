# Imports
import numpy as np
import pyvista as pv
from pathlib import Path
from scipy.io import loadmat

# Variables
INPUT_PATH   = Path(r"data\FemurInput.mat")
RESULTS_PATH = Path(r"results")

SCENARIOS = {
    "Sideways Fall":     RESULTS_PATH / "sideways_fall_results.mat",
    "Single Leg Stance": RESULTS_PATH / "single_leg_stance_results.mat",
}


# Functions
def inspect_mesh(mesh, bc_coords, title_str):
    pl = pv.Plotter()
    pl.set_background("white")
    pl.add_mesh(mesh, color=[0.4, 0.8, 0.7], opacity=0.6,
                show_edges=True, edge_color="gray")
    pl.add_points(pv.PolyData(bc_coords), color="blue",
                  point_size=8, render_points_as_spheres=True)
    pl.add_text(f'Input mesh – "{title_str}"', font_size=11, color="black")
    pl.add_axes(xlabel="X [mm]", ylabel="Y [mm]", zlabel="Z [mm]")
    pl.show()


def inspect_emod(mesh, E_vals, title_str):
    m = mesh.copy()
    m.cell_data["E [MPa]"] = E_vals
    pl = pv.Plotter()
    pl.set_background("white")
    pl.add_mesh(m, scalars="E [MPa]", cmap="jet",
                scalar_bar_args={"title": "E [MPa]", "color": "black"})
    pl.add_text(f'Young\'s modulus – "{title_str}"', font_size=11, color="black")
    pl.add_axes(xlabel="X [mm]", ylabel="Y [mm]", zlabel="Z [mm]")
    pl.show()


def inspect_outputs(mesh, def_mag_elem, epsi1, epsi3, title_str):
    pl = pv.Plotter(shape=(1, 3), window_size=(2800, 850))
    pl.set_background("white")

    datasets = [
        (def_mag_elem, "||u|| [mm]",  "Deformation magnitude"),
        (epsi1,        "e1",          "Max principal strain e1"),
        (epsi3,        "e3",          "Min principal strain e3"),
    ]

    for col, (scalars, bar_title, subplot_title) in enumerate(datasets):
        pl.subplot(0, col)
        m = mesh.copy()
        m.cell_data[bar_title] = scalars
        pl.add_mesh(m, scalars=bar_title, cmap="jet",
                    scalar_bar_args={"title": bar_title, "color": "black"})
        pl.add_text(subplot_title, font_size=10, color="black")
        pl.add_axes(xlabel="X [mm]", ylabel="Y [mm]", zlabel="Z [mm]")

    pl.add_text(f'Outputs – "{title_str}"', font_size=13,
                color="black", position="upper_edge")
    pl.show()


# Main
inp      = loadmat(INPUT_PATH)
nodes    = inp["nodes"]
elements = inp["elements"]

coords  = nodes[:, 1:4].astype(float)
conn    = (elements[:, 1:5] - 1).astype(int)
E_vals  = elements[:, 5].astype(float)
n_nodes = nodes.shape[0]

n_elem     = conn.shape[0]
cells      = np.hstack([np.full((n_elem, 1), 4, dtype=np.int64), conn]).ravel()
cell_types = np.full(n_elem, pv.CellType.TETRA, dtype=np.uint8)
mesh       = pv.UnstructuredGrid(cells, cell_types, coords)

for title_str, mat_path in SCENARIOS.items():
    res   = loadmat(mat_path)
    bcs   = res["bcs"]
    U     = res["U"].ravel()
    epsi1 = res["epsi1"].ravel()
    epsi3 = res["epsi3"].ravel()

    bc_node_ids  = bcs[:, 0].astype(int) - 1
    bc_coords    = coords[bc_node_ids]

    U_mat        = U.reshape(n_nodes, 3)
    def_mag      = np.linalg.norm(U_mat, axis=1)
    def_mag_elem = def_mag[conn].mean(axis=1)

    inspect_mesh(mesh, bc_coords, title_str)
    inspect_emod(mesh, E_vals, title_str)
    inspect_outputs(mesh, def_mag_elem, epsi1, epsi3, title_str)