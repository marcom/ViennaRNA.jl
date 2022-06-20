# Note: ViennaRNA's arrays use 1:n indexing in C. For ViennaRNA's 1:n
#       index range in C we must use a 2:n+1 range when accessing from
#       Julia.
#
# indexfn: calculate index from i,j
# Example for exp_matrices->probs (uses iindx)
#     iindx = fc.uptr.iindx[]
#     indexfn = (i,j) -> iindx[i + 1] - j + 1
# Example for all other matrices which use jindx
#     jindx = fc.uptr.jindx[]
#     indexfn = (i,j) -> jindx[j + 1] + i + 1)

function unsafe_loadmat(fc::FoldCompound, ptr::UnsafePtr{T}; indexfn::Function) where {T}
    ptr == C_NULL && return nothing
    len = Int(fc.uptr.length[])
    @assert len == length(fc)
    m = zeros(T, len, len)
    for i = 1:len, j = i:len
        m[i, j] = ptr[indexfn(i,j)]
    end
    return m
end

function unsafe_loadvec(fc::FoldCompound, ptr::UnsafePtr{T}) where {T}
    ptr == C_NULL && return nothing
    len = Int(fc.uptr.length[])
    @assert len == length(fc)
    # note +1 in array indexing
    v = [T(ptr[i + 1]) for i = 1:len]
    return v
end
