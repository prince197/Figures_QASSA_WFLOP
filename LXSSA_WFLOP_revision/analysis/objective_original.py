"""Authors' objective.py, supplied September 2026 (kept verbatim apart from this docstring and
condensed whitespace). Used as the reference that authors_objective.py is checked against."""
import numpy as np

R = 38.5
K = 0.075
CT = 0.8

A = 1 - np.sqrt(1 - CT)
B = K / R
ALPHA = np.arctan(K)

PENALTY = 1e10

IDEAL_POWER_SCEN1 = 14045.7374

OMEGA_SCEN1 = np.array([
    0, 0.01, 0.01, 0.01, 0.01,
    0.20,
    0.60,
    0.01, 0.01, 0.01, 0.01,
    0.01, 0.01, 0.01, 0.01,
    0.01, 0.01, 0.01, 0.01,
    0.01, 0.01, 0.01, 0.01,
    0
])

OMEGA_SCEN2 = np.array([
    0.0002, 0.0080, 0.0227, 0.0242, 0.0225, 0.0339, 0.0423, 0.0290,
    0.0617, 0.0813, 0.0994, 0.1394, 0.1839, 0.1115, 0.0765, 0.0080,
    0.0051, 0.0019, 0.0012, 0.0010, 0.0017, 0.0031, 0.0097, 0.0317
])

PSI_SCEN2 = np.array([
    7.0, 5.0, 5.0, 5.0, 5.0, 4.0, 5.0, 6.0, 7.0, 7.0, 8.0, 9.5,
    10.0, 8.5, 8.5, 6.5, 4.6, 2.6, 8.0, 5.0, 6.4, 5.2, 4.5, 3.9
])

PSI_SCEN1 = np.full(24, 13.0)

K_SHAPE = 2

LAMBDA = 140.86
ETTA = -500
P_RATED = 1500

CUT_IN = 3.5
RATED = 14

SPEED = np.array([3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11,
                  11.5, 12, 12.5, 13, 13.5, 14])


def velocity_deficit(positions, wind_dir_deg, base_speed):
    positions = np.asarray(positions, dtype=float)
    theta = np.deg2rad(wind_dir_deg)
    dx = positions[:, 0][:, None] - positions[:, 0][None, :]
    dy = positions[:, 1][:, None] - positions[:, 1][None, :]
    proj = dx * np.cos(theta) + dy * np.sin(theta)
    d = np.abs(proj)
    denom = np.sqrt((dx + (R / K) * np.cos(theta)) ** 2
                    + (dy + (R / K) * np.sin(theta)) ** 2)
    arg = np.clip((proj + R / K) / denom, -1.0, 1.0)
    beta = np.arccos(arg)
    deficit = A / ((1.0 + B * d) ** 2)
    mask = beta < ALPHA
    np.fill_diagonal(mask, False)
    wake_sum = np.sqrt(np.sum((deficit ** 2) * mask, axis=1))
    return base_speed * (1.0 - wake_sum)


def boundary_penalty(positions, farm_radius):
    penalty = 0.0
    N = len(positions)
    for i in range(N):
        g = (positions[i, 0] ** 2 + positions[i, 1] ** 2 - farm_radius ** 2)
        if g > 0:
            penalty += (1 + PENALTY * g) ** 2
    return penalty


def spacing_penalty(positions):
    penalty = 0.0
    N = len(positions)
    min_dist = 8 * R
    for i in range(N):
        for j in range(i + 1, N):
            d = np.linalg.norm(positions[i] - positions[j])
            g = min_dist - d
            if g > 0:
                penalty += (1 + PENALTY * g) ** 2
    return penalty


def expected_power(c, omega=None):
    if omega is None:
        omega = OMEGA_SCEN1
    c = np.asarray(c, dtype=float)
    omega = np.asarray(omega, dtype=float)
    s = SPEED[:, None, None]
    E = np.exp(-(s / c[None, :, :]) ** K_SHAPE)
    diff = E[:-1] - E[1:]
    w = (15.0 * omega)[None, :, None]
    inner = np.sum(w * diff, axis=1)
    mid = ((SPEED[:-1] + SPEED[1:]) / 2.0)[:, None]
    first = np.sum(mid * inner, axis=0)
    Erat = np.exp(-(RATED / c) ** K_SHAPE)
    Ecin = np.exp(-(CUT_IN / c) ** K_SHAPE)
    second = np.sum(15.0 * omega[:, None] * Erat, axis=0)
    third = np.sum(15.0 * omega[:, None] * (Ecin - Erat), axis=0)
    return float(np.sum(LAMBDA * first + P_RATED * second + ETTA * third))


IDEAL_POWER_SCEN2 = float(expected_power(PSI_SCEN2.reshape(24, 1), OMEGA_SCEN2))


def scenario_objective(x, farm_radius, dataset=1):
    if dataset == 2:
        omega, psi, ideal = OMEGA_SCEN2, PSI_SCEN2, IDEAL_POWER_SCEN2
    else:
        omega, psi, ideal = OMEGA_SCEN1, PSI_SCEN1, IDEAL_POWER_SCEN1
    N = len(x) // 2
    positions = x.reshape((N, 2))
    c = []
    for k, direction in enumerate(range(0, 360, 15)):
        c.append(velocity_deficit(positions, direction + 7.5, psi[k]))
    c = np.array(c)
    ep = expected_power(c, omega)
    wake_loss = ideal * N - ep
    wake_loss += boundary_penalty(positions, farm_radius)
    wake_loss += spacing_penalty(positions)
    return wake_loss
