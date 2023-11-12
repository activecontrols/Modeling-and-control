function updateVisualizerData(tout,yout)
    addpath("./visualization_resources/")
    model = "Visualizer";
    load_system(model);
    workspace = get_param(model, "modelworkspace");
    workspace.assignin("sim_output", [tout, yout]);
    save_system(model);
end

