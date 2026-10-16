# PIDC baseline (Chan, Stumpf & Babtie 2017; NetworkInference.jl, as in BEELINE).
# Reads a genes x samples TSV (first column = gene names, header row), writes
# gene1 / gene2 / weight for every pair.
#   ~/.juliaup/bin/julia +1.6 tools/pidc.jl in.tsv out.tsv   (NetworkInference.jl 0.1 resolves only on Julia 1.6)
using NetworkInference
inp, out = ARGS[1], ARGS[2]
nodes = get_nodes(inp; delim = '\t', discretizer = "bayesian_blocks")
net = InferredNetwork(PIDCNetworkInference(), nodes)
write_network_file(out, net)
println("pidc: ", length(nodes), " genes -> ", length(net.edges), " edges")
