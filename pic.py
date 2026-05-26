# ---------------------------------------------------------------------------
#                   SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ---------------------------------------------------------------------------
import math

l0 = 2.*math.pi
t0 = l0
Lsim = [64.*l0, 64.*l0]
Tsim = 1000.*t0
resx = 80.
rest = 177.
# 移除移动窗口相关参数
# time_start_moving_window = Lsim[0]

Main(
    geometry = "2Dcartesian",
    interpolation_order = 2,
    cell_length = [l0/resx, l0/resx],
    grid_length = Lsim,
    number_of_patches = [64, 64],
    timestep = t0/rest,
    number_of_pml_cells = [[10, 10], [100, 100]],
    simulation_time = Tsim,
    EM_boundary_conditions = [
        ["PML", "PML"],
        ["PML", "PML"],
    ],
    solve_poisson = False,
    print_every = 100,
    random_seed = smilei_mpi_rank
)

# 移除移动窗口设置，让模拟窗口固定不动
# MovingWindow(
#     time_start = time_start_moving_window,
#     velocity_x = 1.0
# )

LoadBalancing(
    initial_balance = False,
    every = 50000,
    cell_load = 1.,
    frozen_particle_load = 0.1
)

# 激光参数
box_side = "xmin"
a0 = 1.
omega = 1.
focus = [48.*l0, 16.*l0]  # 调整焦点位置到模拟区域中心
waist = 3.0*l0
incidence_angle = 0./180.*math.pi
polarization_phi = 0.
ellipticity = 0.
time_envelope = tgaussian(fwhm=10.5*t0, center=31.5*t0)
phaseZero = 0.

[dephasing, amplitudeY, amplitudeZ] = transformPolarization(polarization_phi, ellipticity)
amplitudeY = a0*omega*amplitudeY
amplitudeZ = a0*omega*amplitudeZ

Zr = omega*waist**2/2.
Y1 = focus[1]
w = math.sqrt(1./(1.+(focus[0]/Zr)**2))
invwaist2 = (w/waist)**2
coeff = -omega * focus[0]*w**2/(2.*Zr**2)

def spatial(y):
    return math.sqrt(w) * math.exp(-invwaist2*(y-focus[1])**2)

def phase(y):
    return coeff * (y-focus[1])**2

def chirp(t):
    b = 0.
    if t < 200.*t0:
        return 1. - (b/(omega*2))*t + (b*31.5*t0)/omega
    return 1.0

# 等离子体密度参数 - 调整位置到模拟区域中心
Ley0 = -32.*l0  # 调整为更靠近中心
Lr = 1.*l0
Lp = 0.1*l0
n0 = 100
L0 = 1.2*l0

def dense(x,y):
       if (y-Ley0)-(x)<=0 and (y-Ley0)-(x)>-Lr:
             return n0
       elif (y-Ley0)-(x)<=L0 and (y-Ley0)-(x)>0:
             n1 = n0*math.exp(-((y-Ley0)-(x))/Lp)
             return n1
       else:
             return 0.

def feo(x,y):
       if (y-Ley0)-(x)<=0 and (y-Ley0)-(x)>-Lr:
             return n0
       else:
             return 0.


def fet(x,y):
       if (y-Ley0)-(x)<=L0 and (y-Ley0)-(x)>0:
             n1 = n0*math.exp(-((y-Ley0)-(x))/Lp)
             return n1

       else:
             return 0.

# 激光定义
Laser(
    box_side = box_side,
    omega = omega,
    chirp_profile = chirp,
    time_envelope = time_envelope,
    space_envelope = [lambda y: amplitudeZ*spatial(y), lambda y: amplitudeY*spatial(y)],
    phase = [lambda y: phase(y)-phaseZero+dephasing, lambda y: phase(y)-phaseZero],
    delay_phase = [0., dephasing]
)

# 定义整个模拟区域为等离子体诊断区域
plasma_region_xmin = 0.*l0
plasma_region_xmax = 64.*l0
plasma_region_ymin = 0.*l0
plasma_region_ymax = 64.*l0

# 只追踪模拟区域内的粒子
def plasma_region_filter(particles):
    x_condition = (particles.x >= plasma_region_xmin) & (particles.x <= plasma_region_xmax)
    y_condition = (particles.y >= plasma_region_ymin) & (particles.y <= plasma_region_ymax)
    return x_condition & y_condition

Species(
    name = 'ion',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    ionization_model = 'none',
    particles_per_cell = 4,
    c_part_max = 1.0,
    mass = 1836.0,
    charge = 1.0,
    number_density = dense,
    time_frozen = Tsim,
    boundary_conditions = [
        ["remove", "remove"],
        ["remove", "remove"],
    ],
)

Species(
    name = 'eon',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    ionization_model = 'none',
    particles_per_cell = 25,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    number_density = fet,
    time_frozen = 0.,
    boundary_conditions = [
        ["remove", "remove"],
        ["remove", "remove"],
    ],
)

Species(
    name = 'eon1',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    ionization_model = 'none',
    particles_per_cell = 4,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    number_density = feo,
    time_frozen = 0.,
    boundary_conditions = [
        ["remove", "remove"],
        ["remove", "remove"],
    ],
)

Checkpoints(
    dump_step = 20*rest,
    keep_n_dumps = 20,
    exit_after_dump = False,
)

# 诊断整个模拟区域的场
DiagFields(
    every = 300,
    fields = ['Ex','Ey','Bz','Rho_eon'],
)

# 诊断整个模拟区域的电子密度分布
#DiagParticleBinning(
#    deposited_quantity = "weight",
#    every = 20*rest,
#    species = ["eon", "eon1"],
#    axes = [
#        ["x", 0.*l0, 64.*l0, 200],
#        ["y", 0.*l0, 64.*l0, 200]
#    ],
#    time_average = 1
#)

# 诊断整个模拟区域的离子密度分布
#DiagParticleBinning(
#    deposited_quantity = "weight",
#    every = 20*rest,
#    species = ["ion"],
 #   axes = [
#        ["x", 0.*l0, 64.*l0, 100],
#        ["y", 0.*l0, 64.*l0, 100]
#    ],
#    time_average = 1
#)

# 诊断能量标量
#DiagScalar(
#    every = 20*rest,
#    vars = ['Uelm', 'Ukin', 'Utot', 'Ukin_eon', 'Ukin_eon1', 'Ukin_ion']
#)

# 只追踪模拟区域内的电子粒子
DiagTrackParticles(
    species = "eon",
    every = 300,
    flush_every = 2,
    filter = plasma_region_filter,
    attributes = ["x", "y", "px", "py", "pz", "weight"]
)

# 只追踪模拟区域内的电子粒子（第二种）
DiagTrackParticles(
    species = "eon1",
    every = 300,
    flush_every = 2,
    filter = plasma_region_filter,
    attributes = ["x", "y", "px", "py", "pz", "weight"]
)



#DiagPerformances(
#    every = 2*rest,
#)