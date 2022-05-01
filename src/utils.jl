function prob_of_basepairs(bppm::Matrix, pt::Pairtable)
    n = length(pt)
    size(bppm) == (n,n) ||
        throw(ArgumentError("size mismatch for basepair matrix and structure"))
    p_unpaired = 1 .- reshape(sum(bppm; dims=1), (n,))
    return [pt[i] == 0 ? p_unpaired[i] : bppm[i, pt[i]] for i = 1:n]
end

prob_of_basepairs(bppm::Matrix, structure::AbstractString) =
    prob_of_basepairs(bppm, Pairtable(structure))

prob_of_basepairs(fc::FoldCompound, pt::Pairtable) =
    prob_of_basepairs(bpp(fc), pt)

prob_of_basepairs(fc::FoldCompound, structure::AbstractString) =
    prob_of_basepairs(bpp(fc), Pairtable(structure))

prob_of_basepairs(sequence::AbstractString, pt::Pairtable) =
    prob_of_basepairs(bpp(sequence), pt)

prob_of_basepairs(sequence::AbstractString, structure::AbstractString) =
    prob_of_basepairs(bpp(sequence), Pairtable(structure))
