from repliclade.sim.simulator import Simulator


def __prompt_rerun():
    rerun = input("Would you like to rerun Repliclade again? Y or N: ").lower()
    while rerun not in ["y", "n"]:
        rerun = input(
            "Invalid Input. Would you like to rerun Repliclade again? Y or N: "
        ).lower()

    return True if rerun == "y" else False


def main():
    simulator = Simulator()
    simulator.run_simulation()
    rerun = __prompt_rerun()
    while rerun is True:
        simulator = Simulator()
        simulator.run_simulation()
        rerun = __prompt_rerun()
