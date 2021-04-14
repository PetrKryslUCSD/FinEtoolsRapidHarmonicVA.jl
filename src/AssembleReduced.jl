#__precompile__(true)

using LinearAlgebra
using FinEtools
using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.AssemblyModule: AbstractSysmatAssembler


mutable struct SysmatAssemblerReduced{T<:Number} <: AbstractSysmatAssembler
    m::Matrix{T} # reduced system matrix
    ndofs_row::FInt; ndofs_col::FInt;
    red_ndofs_row::FInt; red_ndofs_col::FInt;
    t::Matrix{T} # transformation matrix
end

function SysmatAssemblerReduced(t::Matrix{T}) where {T<:Number}
    ndofs_row = ndofs_col = size(t, 1)
    red_ndofs_row = red_ndofs_col = size(t, 2)
    m = fill(zero(T), red_ndofs_row, red_ndofs_col)
    return SysmatAssemblerReduced{T}(m, ndofs_row, ndofs_col, red_ndofs_row, red_ndofs_col, t)
end


function startassembly!(self::SysmatAssemblerReduced{T}, elem_mat_nrows::FInt, elem_mat_ncols::FInt, elem_mat_nmatrices::FInt, ndofs_row::FInt, ndofs_col::FInt) where {T<:Number}
    @assert self.ndofs_row == ndofs_row
    @assert self.ndofs_col == ndofs_col
    self.m[:] .= zero(T)
    return self
end

function assemble!(self::SysmatAssemblerReduced{T}, mat::FMat{T}, dofnums_row::FIntVec, dofnums_col::FIntVec) where {T<:Number}
    R, C = size(mat)
    @assert R == C "The matrix must be square"
    @assert dofnums_row == dofnums_col "The degree of freedom numbers must be the same for rows and columns"
    for i in 1:R
        gi = dofnums_row[i]
        if gi < 1
            mat[i, :] = zero(T)
            dofnums_row[i] = 1
        end
    end
    for i in 1:C
        gi = dofnums_col[i]
        if gi < 1
            mat[:, i] = zero(T)
            dofnums_col[i] = 1
        end
    end
    lt = self.t[]
    self.m .+= lt' * mat * lt
    return self
end

function assemble!(self::SysmatAssemblerReduced{T}, mat::FMat{T}, dofnums_row::FIntMat, dofnums_col::FIntMat) where {T<:Number}
    return assemble!(self, mat, vec(dofnums_row), vec(dofnums_col))
end

"""
    makematrix!(self::SysmatAssemblerReduced)

Make a sparse matrix.
"""
function makematrix!(self::SysmatAssemblerReduced)
    return self.m
end


