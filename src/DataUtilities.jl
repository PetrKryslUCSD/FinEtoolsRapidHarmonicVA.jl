using JSON
using DelimitedFiles
using SparseArrays
using HDF5
using FinEtools

const _sep = "-"

function getcleanname(name)
    cleanname = replace(replace(replace(name,'.'=>'_'),':'=>'_'),' '=>'_')    
    return cleanname
end

function with_extension(filename, ext)
    if match(Regex(".*\\." * ext * "\$"), filename) == nothing
        filename = filename * "." * ext
    end
    return filename
end

function retrieve_json(j)
    j = with_extension(j, "json")
    return open(j, "r") do file
        JSON.parse(file)
    end
end

function store_json(j, d)
    j = with_extension(j, "json")
	open(j, "w") do file
		JSON.print(file, d, 4)
	end
end

const MATRIX_UNKNOWN = -1
const MATRIX_DENSE = 0
const MATRIX_SPARSE = 1

"""
    retrieve_matrix(fname)

Retrieve a matrix.
"""
function retrieve_matrix(fname)
    fname = with_extension(fname, "h5")
    typ = MATRIX_UNKNOWN
    f = try
        h5open(fname, "r") 
    catch SystemError
        nothing
    end
    typ = read(f, "matrix_type")
    if typ == MATRIX_DENSE
        return read(f, "matrix")
    elseif typ == MATRIX_SPARSE
        I = read(f, "I")
        J = read(f, "J")
        V = read(f, "V")
        nrows = read(f, "nrows")
        ncols = read(f, "ncols")
        return sparse(I, J, V, nrows, ncols)
    end
    return nothing
end

"""
    store_matrix(fname, matrix)

Store a matrix.
"""
function store_matrix(fname, matrix)
    fname = with_extension(fname, "h5")
    h5open(fname, "w") do file
        write(file, "matrix_type", MATRIX_DENSE) 
        write(file, "matrix", matrix) 
    end
end


"""
    store_matrix(fname, matrix)

Store a matrix.
"""
function store_matrix(fname, matrix::SparseArrays.SparseMatrixCSC{Float64, Int64})
    I, J, V = findnz(matrix)
    fname = with_extension(fname, "h5")
    h5open(fname, "w") do file
        write(file, "matrix_type", MATRIX_SPARSE) 
        write(file, "I", I) 
        write(file, "J", J) 
        write(file, "V", V) 
        write(file, "nrows", size(matrix, 1)) 
        write(file, "ncols", size(matrix, 2)) 
    end
end



