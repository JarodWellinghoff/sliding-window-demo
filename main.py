import os
import pydicom
import numpy as np
from glob import glob
from scipy.optimize import minimize
import matplotlib.pyplot as plt

def load_ordered_slices(series_dir: str):
    paths = glob(os.path.join(series_dir, "*"))
    ds_list = [pydicom.dcmread(p, stop_before_pixels=True) for p in paths]

    # keep only images from the same series and not localizers
    siu0 = ds_list[0].SeriesInstanceUID
    imgs = [ds for ds in ds_list
            if getattr(ds, "SeriesInstanceUID", None) == siu0 and
               "ImageType" in ds and "LOCALIZER" not in [t.upper() for t in ds.ImageType]]

    if not imgs:
        imgs = ds_list  # fallback if ImageType missing

    # try geometric sort via IOP + IPP
    try:
        io = np.array([float(x) for x in imgs[0].ImageOrientationPatient])
        r, c = io[:3], io[3:]
        n = np.cross(r, c); n /= np.linalg.norm(n)
        svals = []
        for ds in imgs:
            p = np.array(ds.ImagePositionPatient, dtype=float)
            s = float(np.dot(p, n))
            svals.append(s)
        order = np.argsort(svals)

    except Exception:
        # Fallback: InstanceNumber, then AcquisitionTime, then SOPInstanceUID
        def sort_key(ds):
            inn = getattr(ds, "InstanceNumber", None)
            at  = getattr(ds, "AcquisitionTime", "")
            uid = getattr(ds, "SOPInstanceUID", "")
            return (0, inn) if isinstance(inn, int) else (1, at or uid)
        imgs.sort(key=sort_key)
        return imgs  # best-effort fallback

    return [imgs[i] for i in order]

def slice_positions(series_dir: str):
    ds_list = load_ordered_slices(series_dir)

    # Orientation → normal
    io = np.array([float(x) for x in ds_list[0].ImageOrientationPatient])
    r, c = io[:3], io[3:]
    n = np.cross(r, c); n = n / np.linalg.norm(n)

    # Centers (3D) and scalar positions along n
    centers = np.array([np.array(ds.ImagePositionPatient, float) for ds in ds_list])
    scalars = centers @ n  # mm

    # Sort geometrically
    order = np.argsort(scalars)
    scalars = scalars[order]
    ds_sorted = [ds_list[i] for i in order]

    # Thickness (fallback to median)
    thk = np.array([float(getattr(ds, "SliceThickness", np.nan)) for ds in ds_sorted])
    if np.isnan(thk).all():
        raise ValueError("SliceThickness missing in all slices.")
    thk = np.where(np.isnan(thk), np.nanmedian(thk), thk)

    # Scalar positions relative to first slice center
    s0 = scalars[0]
    scalars_rel = scalars - s0  # mm

    # Pack results keyed by InstanceNumber (fallback to SOP)
    out = {}
    for ds, s_rel, t in zip(ds_sorted, scalars_rel, thk,):
        key = getattr(ds, "InstanceNumber", None)
        if key is None:
            key = f"SOP:{ds.SOPInstanceUID}"
        out[key] = {
            "rel_center_mm": float(s_rel),        # along-stack coordinate (0 at first slice center)
            "thickness_mm": float(t),
        }

    return out

def coverage_mm(start, end=None):
    # Compute the coverage between two slices
    if end is None:
        return start["thickness_mm"]
    distance = distance_mm(start, end)
    thickness = thickness_mm(start, end)
    return distance + thickness

def distance_mm(start, end=None):
    # Compute the distance between two slices
    if end is None:
        return 0
    start_pos = start["rel_center_mm"]
    end_pos = end["rel_center_mm"]
    return np.abs(end_pos - start_pos)

def thickness_mm(start, end=None):
    # Compute the thickness between two slices
    if end is None:
        return start["thickness_mm"]
    return (start["thickness_mm"] + end["thickness_mm"]) / 2

def get_future_positions(positions, start_index):
    future_positions_mm = {
        key: pos for key, pos in positions.items() if key > start_index
    }
    return future_positions_mm

def get_asdf(positions, start_index, target_mm):
    future_positions_mm = get_future_positions(positions, start_index)
    coverages_mm = {
        key: coverage_mm(positions[start_index], pos) for key, pos in future_positions_mm.items()
    }
    errors = {key: (abs(c - target_mm) / target_mm) * 100 for key, c in coverages_mm.items()}
    min_error_key = min(errors.items(), key=lambda x: x[1])[0]
    min_error = errors[min_error_key]
    min_error_coverage = coverages_mm[min_error_key]
    return min_error_key, min_error, min_error_coverage

def get_asdf1(positions, start_index, target_mm):
    future_positions_mm = get_future_positions(positions, start_index)
    distances_mm = {
        key: distance_mm(positions[start_index], pos) for key, pos in future_positions_mm.items()
    }
    errors = {key: (abs(c - target_mm) / target_mm) * 100 for key, c in distances_mm.items()}
    min_error_key = min(errors.items(), key=lambda x: x[1])[0]
    min_error = errors[min_error_key]
    min_error_distance = distances_mm[min_error_key]
    return min_error_key, min_error, min_error_distance

if __name__ == "__main__":
    import json, os
    series_dir = os.path.join("W:", "L067_FD_1_0_B30F_0001")
    positions = slice_positions(series_dir)
    smallest_key = min(positions.keys())
    largest_key = max(positions.keys())

    targets = {
        "total_coverage": coverage_mm(positions[smallest_key], positions[largest_key]),
        "window_coverage": 150,
        "step_size": 60
    }

    raw_weights = {
        "total_coverage": 0.5,
        "window_coverage": 0.3,
        "step_size": 0.2
    }
    total_weight = sum(raw_weights.values())
    if total_weight == 0:
        weight_count = len(raw_weights)
        weights = {key: 1 / weight_count for key in raw_weights}
    else:
        weights = {key: val / total_weight for key, val in raw_weights.items()}

    pt = [targets["window_coverage"], targets["step_size"]]

    def get_windows(x):
        window_start_img = smallest_key
        window_end_img, e1, c = get_asdf(positions, window_start_img, x[0])
        next_window_start_img, e2, d = get_asdf1(positions, window_start_img, x[1])
        windows = [{
                "window": {
                    "start_img": window_start_img,
                    "start_idx": window_start_img - 1,
                    "end_img": window_end_img,
                    "end_idx": window_end_img - 1,
                    "error": e1,
                    "coverage_mm": c,
                    "coverage_imgs": window_end_img - window_start_img + 1
                },
                "step": {
                    "img": next_window_start_img,
                    "idx": next_window_start_img - 1,
                    "error": e2,
                    "distance_mm": d,
                    "distance_imgs": next_window_start_img - window_start_img + 1
                }
            }]
        while True:
            window_start_img = next_window_start_img
            window_end_img, e1, c = get_asdf(positions, window_start_img, x[0])
            if e1 > 1:
                break
            next_window_start_img, e2, d = get_asdf1(positions, window_start_img, x[1])
            window = {
                "window": {
                    "start_img": window_start_img,
                    "start_idx": window_start_img - 1,
                    "end_img": window_end_img,
                    "end_idx": window_end_img - 1,
                    "error": e1,
                    "coverage_mm": c,
                    "coverage_imgs": window_end_img - window_start_img + 1
                },
                "step": {
                    "img": next_window_start_img,
                    "idx": next_window_start_img - 1,
                    "error": e2,
                    "distance_mm": d,
                    "distance_imgs": next_window_start_img - window_start_img + 1
                }
            }
            windows.append(window)
        return windows

    def get_stats(windows):
        first_window = windows[0]
        last_window = windows[-1]
        coverage = coverage_mm(positions[first_window["window"]["start_img"]], positions[last_window["window"]["end_img"]])
        window_coverage = sum(window["window"]["coverage_mm"] for window in windows) / len(windows)
        step_size = sum(window["step"]["distance_mm"] for window in windows) / len(windows)
        return {
            "total_coverage": coverage,
            "window_coverage": window_coverage,
            "step_size": step_size
        }

    def get_error(stats, key):
        return (abs(stats[key] - targets[key]) / targets[key]) * 100

    def get_errors(stats):
        return {
            key: get_error(stats, key) for key in stats
        }

    def get_cost(errors, key):
        return weights[key] * errors[key]

    def get_cost_sum(errors):
        return sum(get_cost(errors, key) for key in errors)

    def get_costs(errors):
        c = {
            key: get_cost(errors, key) for key in errors
        }
        c["sum"] = get_cost_sum(errors)
        return c

    def get_results(x):
        windows = get_windows(x)
        stats = get_stats(windows)
        errors = get_errors(stats)
        cost = get_costs(errors)
        return {
            "windows": windows,
            "stats": stats,
            "errors": errors,
            "cost": cost
        }

    # objective function
    def objective(x):
        results = get_results(x)
        return results["cost"]["sum"]
    
    def safe_objective(x):
        try:
            val = objective(np.asarray(x, float))
            if not np.isfinite(val):
                return 1e12
            return val
        except Exception:
            return 1e12

    iter_x = []          # parameter vectors per iteration
    iter_cost = []       # objective values per iteration
    iter_xi = []
    iter_costi = []

    def nm_callback(xk):
        iter_xi.append(np.array(xk, dtype=float))
        print(f"Iteration {len(iter_xi)}")
        # guard: objective may throw if xk is infeasible
        try:
            iter_costi.append(float(objective(np.asarray(xk, float))))
        except Exception:
            iter_costi.append(np.inf)

    iterations = 10

    x0 = np.array([targets["window_coverage"], targets["step_size"]], float)
    res = minimize(
            safe_objective, x0, method="Nelder-Mead",
            callback=nm_callback,
            options={
                "disp": True,
                "maxiter": 2000,
                "xatol": 1e-4,
                "fatol": 1e-4,
                "initial_simplex": np.array([x0 * [0.9, 0.9], x0 * [1.1, 0.9], x0 * [0.9, 1.1]])   # tweak if you like
            }
        )
    iter_x.append(iter_xi)
    iter_cost.append(iter_costi)
    iter_xi = []
    iter_costi = []

    for i in range(iterations-1):
        xi = x0 * np.random.uniform(0.5, 1.5, size=(2,))
        res0 = minimize(
            safe_objective, xi, method="Nelder-Mead",
            callback=nm_callback,
            options={
                "disp": True,
                "maxiter": 2000,
                "xatol": 1e-4,
                "fatol": 1e-4,
                "initial_simplex": np.array([xi * [0.9, 0.9], xi * [1.1, 0.9], xi * [0.9, 1.1]]) # tweak if you like
            }
        )
        iter_x.append(iter_xi)
        iter_cost.append(iter_costi)
        iter_xi = []
        iter_costi = []
        if res0.fun < res.fun:
            res = res0


    paths = []
    for x in iter_x:
        paths.append(np.vstack(x))

    # === 1) Convergence plot ===
    plt.figure()
    for cost in iter_cost:
        plt.semilogy(cost)
    plt.xlabel("Iteration")
    plt.ylabel("Cost f(x)")
    plt.title("Nelder–Mead: Objective Convergence")
    plt.tight_layout()

    # === 2) Parameter traces ===
    plt.figure()
    for path in paths:
        plt.plot(path[:,0], label="window_coverage (mm)")
        plt.plot(path[:,1], label="step_size (mm)")
    plt.xlabel("Iteration")
    plt.ylabel("Parameter value")
    plt.title("Parameter Traces")
    # plt.legend()
    plt.tight_layout()

    # === 3) Path on cost contours (since dim=2) ===
    # Choose a sane box around your start; widen if needed
    w0, s0 = np.mean([np.mean(x, axis=0) for x in iter_x], axis=0)
    w_min, s_min = np.min([np.min(x, axis=0) for x in iter_x], axis=0) * 0.9
    w_max, s_max = np.max([np.max(x, axis=0) for x in iter_x], axis=0) * 1.1

    grid_size = 100
    w_grid = np.linspace(max(1, w_min), w_max, grid_size)
    s_grid = np.linspace(max(1, s_min), s_max, grid_size)
    W, S = np.meshgrid(w_grid, s_grid)

    # Vectorized-safe evaluation (mesh may hit infeasible pockets)
    def eval_grid(W, S):
        Z = np.empty_like(W, dtype=float)
        it = 0
        total = W.shape[0] * W.shape[1]
        for i in range(W.shape[0]):
            for j in range(W.shape[1]):
                Z[i, j] = safe_objective([W[i, j], S[i, j]])
                it += 1
                if it % 100 == 0:
                    print(f"Evaluated {it}/{total} points {(it/total):.2%}")
        return Z

    Z = eval_grid(W, S)

    plt.figure()
    cs = plt.contourf(W, S, Z)
    plt.colorbar(cs)
    for path in paths:
        plt.plot(path[:,0], path[:,1], marker="o", linewidth=1)
        plt.scatter([path[0][0]], [path[0][1]], marker="x", s=80, label="start")
    plt.scatter([res.x[0]], [res.x[1]], marker="*", s=120, label="best")
    plt.xlabel("window_coverage (mm)")
    plt.ylabel("step_size (mm)")
    plt.title("Nelder–Mead Path on Cost Contours")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # === 4) (Optional) Quick look at the final window plan ===
    final_windows = get_windows(res.x)
    print("Best x:", res.x, " Cost:", res.fun, "  #Windows:", len(final_windows))

    results = get_results(res.x)
    results['targets'] = targets
    results['weights'] = weights
    print(json.dumps(results, indent=2))