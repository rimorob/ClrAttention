# Memory-bounded, threaded PIDC with the SAME definition as NetworkInference.jl 0.1 (tools/pidc.jl), bit-identical weights.
#   ~/.juliaup/bin/julia +1.6 -t 32 tools/pidc_chunked.jl in.tsv out.f64 [log]
# in.tsv: genes x samples, first column gene names, header row (as tools/pidc.jl). out.f64: the PIDC weight of every gene pair,
# upper triangle, row-major (pair (i, j), i < j, at index (i-1)*G - i*(i-1)/2 + (j-i), 1-based), raw little-endian Float64.
#
# What is the package's and what is changed (docs/plans/scoreboard_competitors_plan.md, section 8):
#   package, unchanged functions: get_nodes (Bayesian blocks), get_joint_probabilities, apply_mutual_information_formula,
#     apply_specific_information_formula, the Gamma context fit/cdf and its normal-CLR fallback (apply_clr_context copied verbatim).
#   redundancy: apply_redundancy_formula(p_z, si1, si2, base) = sum(vec(p_z) .* vec(min.(si1, si2))); here the same arithmetic into a
#     preallocated buffer of exactly the same length (same elementwise products, the same sum() on a Vector{Float64}).
#   bookkeeping changed: (1) specific information of all sources about a target t is one contiguous bins(t) x G matrix SI[t] (the package
#     keeps a heap NodePair per ordered pair, 18 GB at 8,745 genes, and the Distributed version copies it to every worker);
#     (2) the puc score of the pair {a, b} is accumulated by one thread over w = 1..G in ascending order, adding the term of the target
#     with the larger index first, which is exactly the order in which the package's serial i<j<k loop adds into that matrix cell, so the
#     floating-point sums are identical to the serial run and there are no races; (3) the Gamma fit of a gene's score vector is done once
#     per gene (the package refits it for every pair; same deterministic fit); (4) the output is the packed upper triangle.
using NetworkInference
using Base.Threads
using Statistics
const NI = NetworkInference
const DI = NetworkInference.Distributions

# ---- kernels as separate functions (type-stable; no boxed captured variables) ----
function node_pairs!(MI, SI, nodes, est, base)
    # rows are handed out in blocks with a full GC between blocks: the package functions allocate small arrays for every pair, and without
    # the collections the garbage (not the data) would dominate the peak memory
    G = length(nodes); blk = 128
    for r0 in 1:blk:G
        rows = r0:min(r0 + blk - 1, G); next = Atomic{Int}(first(rows))
        @threads for _ in 1:nthreads()
            while true
                i = atomic_add!(next, 1); i > last(rows) && break
                for j in i+1:G
                    p, p1, p2 = NI.get_joint_probabilities(nodes[i], nodes[j], est)
                    mi = NI.apply_mutual_information_formula(p, p1, p2, base)
                    si1 = NI.apply_specific_information_formula(p, p1, p2, 1, base)   # source i about target j
                    si2 = NI.apply_specific_information_formula(p, p2, p1, 2, base)   # source j about target i
                    MI[i, j] = mi; MI[j, i] = mi
                    SI[j][:, i] .= si1; SI[i][:, j] .= si2
                end
            end
        end
        GC.gc()
    end
end

@inline function redundancy(pz::Vector{Float64}, M::Matrix{Float64}, lo::Int, hi::Int, buf::Vector{Float64})
    @inbounds for k in 1:length(pz); buf[k] = pz[k] * min(M[k, lo], M[k, hi]); end
    return sum(buf)
end

function puc_pair(a::Int, b::Int, G::Int, mi::Float64, Sa::Matrix{Float64}, Sb::Matrix{Float64}, pa::Vector{Float64}, pb::Vector{Float64},
                  ba::Vector{Float64}, bb::Vector{Float64})
    s = 0.0
    for w in 1:G
        (w == a || w == b) && continue
        # target b, sources a and w (ascending); the pair (a, b) gets the term of source a
        lo = a < w ? a : w; hi = a < w ? w : a
        red = redundancy(pb, Sb, lo, hi, bb); c = (mi - red) / mi; c = isfinite(c) ? c : zero(c); s += c
        # target a, sources b and w (ascending); the pair (a, b) gets the term of source b
        lo = b < w ? b : w; hi = b < w ? w : b
        red = redundancy(pa, Sa, lo, hi, ba); c = (mi - red) / mi; c = isfinite(c) ? c : zero(c); s += c
    end
    return s
end

function puc_scores!(S, MI, SI, PZ, bins, bufs)
    G = length(PZ); next = Atomic{Int}(1)
    @threads for _ in 1:nthreads()
        buf = bufs[threadid()]
        while true
            a = atomic_add!(next, 1); a > G && break
            for b in a+1:G
                s = puc_pair(a, b, G, MI[a, b], SI[a], SI[b], PZ[a], PZ[b], buf[bins[a]], buf[bins[b]])
                S[a, b] = s; S[b, a] = s
            end
        end
    end
end

function context_weights(S, G)
    fits = Vector{Any}(undef, G); mus = zeros(G); vars = zeros(G); nfail = Atomic{Int}(0); next = Atomic{Int}(1)
    @threads for _ in 1:nthreads()
        while true
            i = atomic_add!(next, 1); i > G && break
            v = vcat(S[1:i-1, i], S[i+1:end, i])
            mus[i] = mean(v); vars[i] = var(v)
            fits[i] = try DI.fit(DI.Gamma, v) catch; atomic_add!(nfail, 1); nothing end
        end
    end
    W = Vector{Float64}(undef, G * (G - 1) ÷ 2); nfb = Atomic{Int}(0); next[] = 1
    @threads for _ in 1:nthreads()
        while true
            i = atomic_add!(next, 1); i > G && break
            base_idx = (i - 1) * G - (i * (i - 1)) ÷ 2
            for j in i+1:G
                score = S[i, j]
                w = try
                    (fits[i] === nothing || fits[j] === nothing) && error("no fit")
                    DI.cdf(fits[i], score) + DI.cdf(fits[j], score)
                catch
                    atomic_add!(nfb, 1)
                    di = score - mus[i]; dj = score - mus[j]
                    sqrt((vars[i] == 0 || di < 0 ? 0 : di^2 / vars[i]) + (vars[j] == 0 || dj < 0 ? 0 : dj^2 / vars[j]))
                end
                W[base_idx + (j - i)] = w
            end
        end
    end
    return W, nfail[], nfb[]
end

function main(inp, outp, logp)
    log = open(logp, "w"); say(s) = (println(log, s); flush(log); println(s); flush(stdout))
    t0 = time()
    nodes = get_nodes(inp; delim = '\t', discretizer = "bayesian_blocks")
    G = length(nodes); est = "maximum_likelihood"; base = 2
    bins = [n.number_of_bins for n in nodes]; maxb = maximum(bins)
    say("nodes $(G), bins min/mean/max $(minimum(bins))/$(round(sum(bins)/G, digits=2))/$(maxb), threads $(nthreads()), read $(round(time()-t0, digits=1)) s")
    MI = zeros(Float64, G, G)
    SI = [zeros(Float64, bins[t], G) for t in 1:G]          # SI[t][:, s] = specific information of source s about target t
    PZ = Vector{Float64}[vec(nodes[t].probabilities) for t in 1:G]
    node_pairs!(MI, SI, nodes, est, base)
    say("node pairs done $(round(time()-t0, digits=1)) s; peak RSS GB $(round(Sys.maxrss()/2^30, digits=2))")
    S = zeros(Float64, G, G)
    bufs = [[zeros(Float64, k) for k in 1:maxb] for _ in 1:nthreads()]
    puc_scores!(S, MI, SI, PZ, bins, bufs)
    say("puc scores done $(round(time()-t0, digits=1)) s; peak RSS GB $(round(Sys.maxrss()/2^30, digits=2))")
    SI = nothing; MI = nothing; GC.gc()
    W, nfail, nfb = context_weights(S, G)
    say("gamma fits: $(nfail) genes without a Gamma fit; weights: $(nfb) pairs through the normal-CLR fallback; $(round(time()-t0, digits=1)) s")
    open(outp, "w") do io; write(io, W); end
    say("wrote $(outp): $(length(W)) pair weights; total $(round(time()-t0, digits=1)) s; peak RSS GB $(round(Sys.maxrss()/2^30, digits=2))")
    close(log)
end
main(ARGS[1], ARGS[2], length(ARGS) >= 3 ? ARGS[3] : ARGS[2] * ".log")
