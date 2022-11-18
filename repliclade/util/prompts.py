from repliclade.common.repliclade_constants import EVOL_MODELS


def file_prompt():
    """Prompts file name from the user"""
    filename = input("Please specify the name of the file: ")
    return filename


def prompt_model():
    """Prompts the Evolutionary Model to use from the user"""
    for model in EVOL_MODELS:
        print(model)
    model = input(
        "Please specify the evolutionary model you would like to use from the ones given above: "
    )
    while model.lower() not in EVOL_MODELS:
        model = input(
            "Invalid input.  Please specify the evolutionary model you would like to use from the ones given above: "
        )
    return model.lower()


def cr_prompt():
    """Prompts user whether conserved regions should be considered"""
    decision = input(
        "Would you like to consider conserved regions previously identified?  Please input Y or N: "
    )
    while decision not in ["Y", "y", "N", "n"]:
        decision = input("Invalid Input.  Please specify Y or N")
    return decision


def align_results_prompt():
    """Prompts the user whether an alignment should be performed"""
    decision = input(
        "Would you like to perform an alignment on the results for the sequence/sequences? Please input 'Y' or 'N': "
    )
    while decision not in ["Y", "y", "N", "n"]:
        decision = input("Invalid Input.  Please specify Y or N")
    return decision


def align_single_or_multiple_prompt():
    """Prompts the user whether one or multiple sequences should be aligned"""
    inp = input(
        "Would you like to align all of the sequences and their results or just a specific sequence? Please specify 'single' or multiple': "
    )
    while inp.lower() not in ["single", "multiple"]:
        inp = input("Invalid input.  Please specify single or multiple")
    return inp


# TODO is this needed?
def which_sequence_align_prompt(seq_ids):
    [print(seq_id.lower()) for seq_id in seq_ids]
    seq_ids = [seq.lower() for seq in seq_ids]
    inp = input(
        "Which sequence would you like to align for every iteration throughout the simulation?: "
    )
    while inp.lower() not in seq_ids:
        inp = input("Invalid input.  Please specify one of the sequence Id's: ")
    return (inp, seq_ids.index(inp.lower()))


def prompt_time(coalescence_time):
    print(
        "We have estimated earlier that the total ancestor coalescence time is {}".format(
            coalescence_time
        )
    )
    print(
        "We can provide an amount of generations for you to run the simulation.  Would you like to enter an amount in years for each generation for this specific sequence?"
    )
    yes_or_no = input("Y or N: ")
    while yes_or_no not in ["N", "n", "Y", "y"]:
        yes_or_no = input("Invalid input. Please enter Y or N: ")
    if yes_or_no in ["Y", "y"]:
        years_per_generation = input(
            "Please enter the amount in years for each generation for this input sequence: "
        )
        try:
            years_per_generation = float(years_per_generation)
        except:
            while not isinstance(years_per_generation, float):
                years_per_generation = input(
                    "Invalid input. Please input a valid number: "
                )
                try:
                    years_per_generation = float(years_per_generation)
                except:
                    continue
        generations = int(coalescence_time / years_per_generation)
        print("The simulation will run for {} generations".format(generations))
        return generations
    else:
        generations = input(
            "Please specify the amount of generations you would like to run the simulation for: "
        )
        try:
            generations = int(generations)
        except:
            while not isinstance(generations, int):
                generations = input(
                    "Invalid input.  Please input the amount of generations you would like to run the simulation for as an integer: "
                )
                try:
                    generations = int(generations)
                except:
                    continue
        return generations


def prompt_mutation_rate():
    """
    Prompts the user for mutation rate
    """
    inp = input(
        "Would you like to estimate the effective pop size with our own input? Please enter y or n: "
    )
    while inp.lower() not in ["y", "n"]:
        inp = input("Invalid input.  Please enter y or n: ")
    if inp.lower() == "y":
        mu = 0
        try:
            mu = float(input("Please enter a value for the mutation rate: "))
        except:
            while not isinstance(mu, float):
                mu = input(
                    "Invalid input.  Please enter a valid number for the mutation rate"
                )
                try:
                    mu = float(mu)
                except:
                    continue
        return mu
    else:
        return None
