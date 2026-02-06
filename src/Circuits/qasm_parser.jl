### qasm_parser.jl
##
# This module contains functions for parsing OpenQASM 2.0 into PauliPropagation circuits.
# OpenQASM is a quantum assembly language used by many quantum computing frameworks.
##
###

"""
    parse_qasm(qasm_string::String; nqubits::Union{Int,Nothing}=nothing) -> Vector{Gate}

Parse an OpenQASM 2.0 string into a PauliPropagation circuit.

# Arguments
- `qasm_string::String`: The OpenQASM 2.0 code to parse
- `nqubits::Union{Int,Nothing}=nothing`: Number of qubits. If not provided, inferred from qreg declarations.

# Returns
- `Vector{Gate}`: A circuit (vector of gates) that can be used with PauliPropagation

# Supported Gates
## Single-qubit gates:
- `h q[i]` → Hadamard gate
- `x q[i]` → Pauli X gate
- `y q[i]` → Pauli Y gate
- `z q[i]` → Pauli Z gate
- `s q[i]` → S gate (phase gate)
- `rx(θ) q[i]` → RX rotation
- `ry(θ) q[i]` → RY rotation
- `rz(θ) q[i]` → RZ rotation

## Two-qubit gates:
- `cx q[i],q[j]` (or `cnot q[i],q[j]`) → CNOT gate
- `cz q[i],q[j]` → CZ gate
- `swap q[i],q[j]` → SWAP gate

# Example
```julia
qasm = \"\"\"
OPENQASM 2.0;
qreg q[2];
h q[0];
cx q[0],q[1];
\"\"\"
circuit = parse_qasm(qasm)
```

# Notes
- Qubit indices in OpenQASM are 0-based but converted to 1-based for PauliPropagation
- Measurement operations and classical register operations are ignored
- Include statements are ignored
- Only gate operations are converted to circuit gates
"""
function parse_qasm(qasm_string::String; nqubits::Union{Int,Nothing}=nothing)
    circuit = Gate[]
    num_qubits = nqubits
    
    # Split into lines and process
    lines = split(qasm_string, '\n')
    
    for line in lines
        # Remove comments and trim whitespace
        line = strip(split(line, "//")[1])
        
        # Skip empty lines
        if isempty(line)
            continue
        end
        
        # Skip header, includes, and declarations we don't need
        if startswith(line, "OPENQASM") || 
           startswith(line, "include") || 
           startswith(line, "creg")
            continue
        end
        
        # Parse qreg to get number of qubits if not provided
        if startswith(line, "qreg")
            if num_qubits === nothing
                m = match(r"qreg\s+\w+\[(\d+)\]", line)
                if m !== nothing
                    num_qubits = parse(Int, m.captures[1])
                end
            end
            continue
        end
        
        # Skip measure operations
        if startswith(line, "measure") || startswith(line, "barrier")
            continue
        end
        
        # Remove trailing semicolon
        line = replace(line, ";" => "")
        
        # Parse gate operations
        gate = _parse_gate_line(line)
        if gate !== nothing
            push!(circuit, gate)
        end
    end
    
    return circuit
end


"""
    parse_qasm_file(filename::String; nqubits::Union{Int,Nothing}=nothing) -> Vector{Gate}

Parse an OpenQASM 2.0 file into a PauliPropagation circuit.

# Arguments
- `filename::String`: Path to the OpenQASM file
- `nqubits::Union{Int,Nothing}=nothing`: Number of qubits. If not provided, inferred from qreg declarations.

# Returns
- `Vector{Gate}`: A circuit (vector of gates) that can be used with PauliPropagation
"""
function parse_qasm_file(filename::String; nqubits::Union{Int,Nothing}=nothing)
    qasm_string = read(filename, String)
    return parse_qasm(qasm_string; nqubits=nqubits)
end


"""
    _parse_gate_line(line::String) -> Union{Gate,Nothing}

Parse a single gate line from OpenQASM.
Returns Nothing if the line is not a recognized gate.
"""
function _parse_gate_line(line::String)
    line = strip(line)
    
    # Pattern for single qubit gate: gate q[i]
    # Pattern for single qubit parametrized gate: gate(param) q[i]
    # Pattern for two qubit gate: gate q[i],q[j]
    
    # Try to match parametrized single-qubit gates (rx, ry, rz)
    m = match(r"^(rx|ry|rz)\s*\(\s*([^)]+)\s*\)\s+\w+\[(\d+)\]", line)
    if m !== nothing
        gate_name = m.captures[1]
        param_str = m.captures[2]
        qubit = parse(Int, m.captures[3]) + 1  # Convert to 1-based indexing
        
        # Parse parameter - it could be a number or expression like pi/2
        param = _parse_parameter(param_str)
        
        if gate_name == "rx"
            return PauliRotation(:X, qubit, param)
        elseif gate_name == "ry"
            return PauliRotation(:Y, qubit, param)
        elseif gate_name == "rz"
            return PauliRotation(:Z, qubit, param)
        end
    end
    
    # Try to match single-qubit Clifford gates
    m = match(r"^([hxyzsHXYZS])\s+\w+\[(\d+)\]", line)
    if m !== nothing
        gate_name = lowercase(m.captures[1])
        qubit = parse(Int, m.captures[2]) + 1  # Convert to 1-based indexing
        
        gate_symbol = _qasm_to_clifford_symbol(gate_name)
        if gate_symbol !== nothing
            return CliffordGate(gate_symbol, qubit)
        end
    end
    
    # Try to match two-qubit gates: cx q[i],q[j] or cnot q[i],q[j]
    m = match(r"^(cx|cnot|cz|swap)\s+\w+\[(\d+)\]\s*,\s*\w+\[(\d+)\]", line)
    if m !== nothing
        gate_name = lowercase(m.captures[1])
        qubit1 = parse(Int, m.captures[2]) + 1  # Convert to 1-based indexing
        qubit2 = parse(Int, m.captures[3]) + 1
        
        if gate_name == "cx" || gate_name == "cnot"
            return CliffordGate(:CNOT, [qubit1, qubit2])
        elseif gate_name == "cz"
            return CliffordGate(:CZ, [qubit1, qubit2])
        elseif gate_name == "swap"
            return CliffordGate(:SWAP, [qubit1, qubit2])
        end
    end
    
    # If we get here, we didn't recognize the gate
    # We could either error or silently skip - for now, skip with a warning
    if !isempty(line)
        @warn "Skipping unrecognized gate: $line"
    end
    
    return nothing
end


"""
    _qasm_to_clifford_symbol(gate_name::String) -> Union{Symbol,Nothing}

Map OpenQASM gate names to PauliPropagation Clifford gate symbols.
"""
function _qasm_to_clifford_symbol(gate_name::String)
    gate_map = Dict(
        "h" => :H,
        "x" => :X,
        "y" => :Y,
        "z" => :Z,
        "s" => :S
    )
    
    return get(gate_map, gate_name, nothing)
end


"""
    _parse_parameter(param_str::AbstractString) -> Float64

Parse a parameter string from OpenQASM.
Supports numeric values and expressions involving π (pi).
"""
function _parse_parameter(param_str::AbstractString)
    param_str = strip(String(param_str))
    
    # Replace common patterns
    param_str = replace(param_str, "pi" => "π")
    
    # Try to evaluate as Julia expression
    try
        return eval(Meta.parse(param_str))
    catch
        # If that fails, try to parse as a simple number
        try
            return parse(Float64, param_str)
        catch
            error("Could not parse parameter: $param_str")
        end
    end
end
