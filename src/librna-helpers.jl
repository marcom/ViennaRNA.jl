# Helper functions for LibRNA functionality
#
# These are essentially C code in Julia and could perhaps be included
# in ViennaRNA.

module LibRNA_Helper
function free_subopt_solutions(ptr::Ptr)
    ptr == C_NULL && return
    i = 1
    while true
        sol = unsafe_load(ptr, i)
        sol.structure == C_NULL && break
        Libc.free(sol.structure)
        i += 1
    end
    Libc.free(ptr)
end

function free_structure_list(ptr::Ptr, num::Integer)
    for i = 1:num
        Libc.free(unsafe_load(ptr, i))
    end
    Libc.free(ptr)
end
end # module LibRNA_Helper
