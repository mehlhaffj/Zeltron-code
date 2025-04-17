import numpy as np

c = 299792458e2
e = 4.8032068e-10
me=9.1093897e-28

def print_quantities(qsd):
    max_symbol_len = 0
    max_desc_len = 0
    for tup in qsd:
        quantity, symbol, desc = tup
        max_symbol_len = max(max_symbol_len, len(symbol))
        max_desc_len = max(max_desc_len, len(desc))
        
    for tup in qsd:
        quantity, symbol, desc = tup
        mystr  = desc.ljust(max_desc_len) + ", "
        mystr += symbol.ljust(max_symbol_len) + " = "
        mystr += str(quantity)
        print(mystr)


if __name__ == "__main__":
    # Nominal Larmor radius: rho0 == me c^2 / e b0
    rho0  = 1.0
    # Field strength on z-axis at z=zmin
    # The point (x,y,z) = (0,0,zmin) will hereafter be referred to as the
    # 'nozzle center'.
    # The 2D surface defined by z=zmin and x^2 + y^2 <= rdisk is referred to as
    # the 'nozzle'.
    b0    = (me * c * c) / (e * rho0)
    # Magnetization of the initially stationary background plasma at the nozzle
    # center
    sigbg = 1.5**4  # 2.0**4  # 3.0**4
    # Temperature of the initially stationary background plasma
    thbg  = 1.0 / np.sqrt(sigbg)
    # Magnetization of the nozzle-injected plasma
    siginj = np.sqrt(sigbg)
    # Nozzle-injected multiplicity (with respect to Goldreich-Julian number
    # density)
    kappa = 6.0
    # Number of cells in x-direction
    Nx = 128

    # The light-cylinder radius
    rlc = 2 * kappa * siginj
    # The nozzle radius
    rnozzle = 0.9 * rlc

    # The initial magnetic monopole is placed at (x,y,z)=(0,0,zmin-rmp)
    # It is outside of the simulation domain
    rmp = rlc

    # The number density of the initially uniform/stationary backgroudn plasma
    nbg  = b0**2 / (4 * np.pi * me * c * c * sigbg)
    # The injected number density at the nozzle center
    ninj = b0**2 / (4 * np.pi * me * c * c * siginj)

    # The non-relativistic skin-depth of the injected plasma at the nozzle
    # center
    dinj = np.sqrt(siginj) * rho0
    # The Debye length of the background plasma
    lambg = np.sqrt(sigbg) * rho0 * np.sqrt(thbg)

    # The critical radius where the background cold magnetization, evaluated
    # using the Michel solution, drops below 1
    rc = rmp / rlc * np.sqrt(sigbg) * (1. / (1. + rmp**2 / rnozzle**2)) * rmp

    # Boundaries of the simulation domain
    xmax = Nx # 4 * rc
    xmin = -xmax
    Lx = xmax - xmin

    dx = Lx / Nx
    dt = 0.99 * dx / (c * np.sqrt(3.))

    quantities_symbols_descriptions = [
        (sigbg, "sigma_bg", "Background cold magnetization"),
        (siginj, "sigma_inj", "Nozzle-injected magnetization"),
        (thbg, "theta_bg", "Background temperature"),
        (kappa, "kappa", "Nozzle-center multiplicity"),
        (rlc, "R_LC", "Light cylinder radius"),
        (rho0, "rho_0", "Nominal Larmor radius"),
        (dinj, "d_inj", "Injected plasma non-rel. skin depth"),
        (lambg, "lambda_bg", "Background plasma Debye length"),
        (rnozzle, "R_nozz", "Nozzle radius"),
        (rnozzle / rlc, "R_nozz / R_LC", ""),
        (rlc / dinj, "R_LC / d_inj", ""),
        (ninj, "n_inj", "Nozzle-center injected num. density"),
        (nbg, "n_bg", "Background num. density"),
        (ninj / nbg, "n_inj / n_bg", ""),
        (rmp, "r_mp", "Initial monopole location z=zmin-r_mp"),
        (rc, "R_c", "Collision radius b/w injected and bg plasma"),
        (rc / rlc, "R_c / R_LC", ""),
        (4 * rc / rlc, "4 * R_c / R_LC", ""),
        (4 * rc / rho0, "4 * R_c / rho_0", ""),
        (4 * rc / dinj, "4 * R_c / d_inj", ""),
        (Lx, "L_x", "Simulation length in x-direction"),
        (Nx, "N_x", "Number of cells in x-direction"),
        (dx, "dx", "L_x / N_x"),
        (dx / dinj, "dx / d_inj", ""),
        (c * dt, "c * dt", "CFL timestep"),
        (rlc / (c * dt), "R_LC / (c * dt)", ""),
        (Lx / (c * dt), "L_x / (c * dt)", ""),
        (Lx / rlc, "L_x / R_LC", ""),
    ]

    print_quantities(quantities_symbols_descriptions)
