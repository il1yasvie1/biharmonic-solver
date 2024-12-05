from src.biharmonic import BiharmonicProblem
import hydra
from src.utils import plot_deformation


@hydra.main(config_path="conf", config_name="config", version_base=None)
def run(cfg):
    problem = BiharmonicProblem(cfg)
    y = problem.solve()
    plot_deformation(y, './outputs/figures/solution.png', equal_aspect=True)
    # output_file = VTKFile("./outputs/solution.pvd", project_output=True)
    # output_file.write(y)


if __name__ == "__main__":
    run()
