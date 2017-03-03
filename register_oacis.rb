repo_dir = File.expand_path(File.dirname(__FILE__))

localhost = Host.find_by_name("localhost")

sim_params = {
  name: "EOS_model",
  command: "#{repo_dir}/run.sh",
  support_input_json: false,
  print_version_command: "cd #{repo_dir} && git describe --always",
  parameter_definitions: [
    {key: "model_type", type: "Integer", default: 0, description: "0:standard model, 1: flat degree, 2: flat link weight, 3: flat degree & weight"}, 
    {key: "m", type: "Integer", default: 10, description: "number of interactions"},
    {key: "t", type: "Integer", default: 30000, description: "simulation time"}
  ],
  description: "EOS model",
  executable_on: [ localhost ]
}

sim_name = sim_params[:name]
if Simulator.where(name: sim_name).exists?
  puts "simulator #{sim_name} already exists" 
else
  sim = Simulator.create!(sim_params)
end

