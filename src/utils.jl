# Note
#
# ViennaRNA's arrays use 1:n indexing in C. For ViennaRNA's C 1:n
# indices we must use a 2:n+1 range when accessing from Julia.  This
# is the cause of the +1 increment when accessing ViennaRNA C arrays,
# e.g. jindex[][j + 1].

function unsafe_loadmat(fc::FoldCompound,
                        mat::Union{UnsafePtr{LibRNA.vrna_mx_mfe_s},
                                   UnsafePtr{LibRNA.vrna_mx_pf_s}},
                        field::Symbol)
    #mat == C_NULL && return nothing
    ptr = getproperty(mat, field)
    return unsafe_loadmat(fc, ptr)
end

function unsafe_loadmat(fc::FoldCompound, ptr::UnsafePtr{Ptr{T}}) where {T}
    ptr[] == C_NULL && return nothing
    len = Int(fc.uptr.length[])
    @assert len == length(fc)
    jindx = fc.uptr.jindx
    m = zeros(T, len, len)
    for i = 1:len
        for j = i:len
            # note +1 in array indexing, C: jindx[j] + i
            ij = jindx[][j + 1] + i
            m[i, j] = ptr[][ij + 1]
        end
    end
    return m
end

function unsafe_loadvec(fc::FoldCompound, mat::UnsafePtr, field::Symbol)
    mat == C_NULL && return nothing
    ptr = getproperty(mat, field)
    return unsafe_loadvec(fc, ptr)
end

function unsafe_loadvec(fc::FoldCompound, ptr::UnsafePtr{Ptr{T}}) where {T}
    ptr[] == C_NULL && return nothing
    len = Int(fc.uptr.length[])
    @assert len == length(fc)
    # note +1 in array indexing
    v = [T(ptr[][i + 1]) for i = 1:len]
    return v
end
