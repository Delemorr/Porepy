import numpy as np
import porepy as pp
from porepy.applications.md_grids.domains import nd_cube_domain
from src.equivalent_permeability_fracture import ekv
from porepy.models.constitutive_laws import DimensionDependentPermeability

from src.fracture_corner import corner_fractur
import matplotlib.pyplot as plt

import pandas as pd
import shapely

aperture = 0.0005
k1 = 100 * 10 ** (-15)
ktr = aperture**2/12
P1 = 10 * 10 ** 4
P2 = 5 * 10 ** 4
m = 0.25

"b - размер стороны трещины"
b = 2
"z-размер ячейки"
z = 1
"ro - плотность"
ro = 1000
"mu - вязкость"
mu = 0.01
ntr = 10
delfi = 3
fi = 0
N = int(90 / delfi)
"L - Длина области"
L = int(b*ntr+b)
print(L)
# ---------------------------------------------------------------------------------------------
"Однако класс SolidConstants по умолчанию не охватывает отдельные трещины и проницаемость матрицы"

class FractureSolidConstants(pp.SolidConstants):
    """Solid constants tailored to the current model."""

    @property
    def default_constants(self):
        """Add the additional parameter `fracture_permeability`."""
        constants = super().default_constants
        constants.update({"fracture_permeability": 1.0})
        return constants

    def fracture_permeability(self) -> float:
        """Permeability of fractures [m^2]."""
        return self.convert_units(self.constants["fracture_permeability"], "m^2")

K = [] # Список экв проницаемости в зависимости от угла
FI = []
for i in range(N + 1):
    # fi = 60
    # ---------------------------------------------------------------------------------------------
    # "Построение системы трещин"
    # box = {'xmin': 0, 'xmax': L, 'ymin': 0, 'ymax': L}
    # domain = pp.Domain(bounding_box=box)
    # fractures = corner_fractur(b, fi, L, ntr)
    # network_2d = pp.create_fracture_network(fractures, domain)
    # mesh_args: dict[str, float] = {"cell_size": 1.0, "cell_size_fracture": 0.9}
    # mdg = pp.create_mdg("simplex", mesh_args, network_2d)
    # "Построение графика"
    # network_2d.plot()
    # ---------------------------------------------------------------------------------------------
    """Класс для изменения геометрии"""
    class Geometry:
        def set_domain(self) -> None:
            """ Определение двумерной квадратной области"""
            size = self.solid.convert_units(L, "m")
            self._domain = nd_cube_domain(2, size)

        def set_fractures(self) -> None:
            """Установка Трещин"""
            frac=corner_fractur(b,fi,L,ntr)
            self._fractures = frac

        def grid_type(self) -> str:
            """ Выбор типа сетки для нашей расчетной области."""
            return self.params.get("grid_type", "simplex")

        def meshing_arguments(self) -> dict:
            """Определение размера ячейки"""
            cell_size = self.solid.convert_units(z, "m")
            mesh_args: dict[str, float] = {"cell_size": cell_size}
            return mesh_args

    # ---------------------------------------------------------------------------------------------
    "Класс для изменения ГУ"
    class BoundaryConditions:

        def bc_type_darcy_flux(self, sd: pp.Grid) -> pp.BoundaryCondition:
            """ Назначение Дирихле западной и восточной границам. Остальные по умолчанию — Нейман."""
            bounds = self.domain_boundary_sides(sd)
            bc = pp.BoundaryCondition(sd,  bounds.west + bounds.east, "dir")
            return bc

        def bc_values_pressure(self, boundary_grid: pp.BoundaryGrid) -> np.ndarray:
            """ Назначение давления на границах."""
            bounds = self.domain_boundary_sides(boundary_grid)
            values = np.zeros(boundary_grid.num_cells)
            values[bounds.west] = self.fluid.convert_units(P1, "Pa")
            values[bounds.east] = self.fluid.convert_units(P2, "Pa")
            return values

    class Permeability(DimensionDependentPermeability):
        params: dict
        solid: FractureSolidConstants

        def fracture_permeability(self, subdomains: list[pp.Grid]) -> pp.ad.Operator:
            size = sum([sd.num_cells for sd in subdomains])
            permeability = pp.wrap_as_dense_ad_array(
                self.solid.fracture_permeability(), size, name="fracture permeability"
            )
            return self.isotropic_second_order_tensor(subdomains, permeability)
    # ---------------------------------------------------------------------------------------------
    """Объединение модели"""
    class Flow_Model(
        Geometry,
        Permeability,
        BoundaryConditions,
        pp.fluid_mass_balance.SinglePhaseFlow,):
        ...
    # ---------------------------------------------------------------------------------------------
    """параметры задачи"""
    solid_constants = FractureSolidConstants(
        {
            "residual_aperture": aperture,
            "fracture_permeability": ktr,
            "porosity": m,
            "permeability": k1
        })
    fluid_constants = pp.FluidConstants({"viscosity": mu, "density": ro})
    params = {"material_constants": {"fluid": fluid_constants, "solid": solid_constants}}

    # ---------------------------------------------------------------------------------------------
    """Запуск модели"""
    model = Flow_Model(params)
    pp.run_time_dependent_model(model, params)
    # pp.plot_grid(model.mdg, "pressure", figsize=(10, 8), linewidth=0.25, title="Pressure distribution", plot_2d=True, pointsize=20, fracturewidth_1d=3)

    # ---------------------------------------------------------------------------------------------
    """Информация про проницаемость"""
    sd = model.mdg.subdomains()
    data = model.mdg.subdomain_data(sd[0])
    gg = data['parameters']['flow']['second_order_tensor'].values
    # pp.plot_grid(sd[0], gg[0,0], title="Проницаемость", plot_2d=True)

    # ---------------------------------------------------------------------------------------------
    "Запись геометрии"
    # g = model.mdg.compute_geometry()
    # f = shapely.to_wkt(g)
    # file = open("geometry"+str(i)+".txt", "w")
    # file.write(str(f))
    # file.close()
    # ---------------------------------------------------------------------------------------------
    """Поиск эквивалентной проницаемости вариант 1"""

    k = ekv(L, ro, mu, P1, P2, model)
    K.append(k)
    FI.append(fi)

    # ---------------------------------------------------------------------------------------------
    "Заполнение тбл результатов"
    data = {
        'Угол': [str(fi)],
        'Эквивалентная проницаемость': [str(k)],
    }
    df = pd.DataFrame(data)
    if fi == 0:
        df.to_csv('ekv.csv', index=False, mode='a', sep='\t' )
    if fi > 0:
        df.to_csv('ekv.csv', index=False, mode='a', sep='\t', header=False )
    fi = fi + delfi

plt.plot(FI,K)
plt.xlabel('fi, угол поворота')
plt.ylabel('k, м^2')
plt.semilogy()
plt.show()

