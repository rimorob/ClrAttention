# Parallel PIDC (NetworkInference.jl @distributed over Distributed workers): `julia +1.6 -p N tools/pidc_par.jl in.tsv out.tsv`. Same algorithm as tools/pidc.jl; edge weights agree with the single-process run to about 1e-6 (116 of 1.35M net1 pairs differ by more than 1e-5, max 0.002; top-100k overlap 1.0).
using Distributed
@everywhere using NetworkInference
inp, out = ARGS[1], ARGS[2]
nodes = get_nodes(inp; delim = '\t', discretizer = "bayesian_blocks")
net = InferredNetwork(PIDCNetworkInference(), nodes)
write_network_file(out, net)
println("pidc: ", length(nodes), " genes -> ", length(net.edges), " edges; workers ", nworkers())
