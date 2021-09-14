def prompt_model():
    evol_models = ["kimura", "jukescantor", "felsenstein", "hasegawa"]
    for model in evol_models:
        print(model)
    model = input(
        "Please specify the evolutionary model you would like to use from the ones given above: "
    )
    while model.lower() not in evol_models:
        model = input(
            "Invalid input.  Please specify the evolutionary model you would like to use from the ones given above: "
        )
    return model.lower()