using Test
using PauliPropagation

@testset "Test OpenQASM Parser" begin
    
    @testset "Basic Single-Qubit Gates" begin
        qasm = """
        OPENQASM 2.0;
        qreg q[2];
        h q[0];
        x q[1];
        """
        circuit = parse_qasm(qasm)
        
        @test length(circuit) == 2
        @test circuit[1] isa CliffordGate
        @test circuit[1].symbol == :H
        @test circuit[1].qinds == [1]  # 0-based -> 1-based
        @test circuit[2] isa CliffordGate
        @test circuit[2].symbol == :X
        @test circuit[2].qinds == [2]
    end
    
    @testset "Pauli Gates" begin
        qasm = """
        OPENQASM 2.0;
        qreg q[3];
        x q[0];
        y q[1];
        z q[2];
        """
        circuit = parse_qasm(qasm)
        
        @test length(circuit) == 3
        @test circuit[1].symbol == :X
        @test circuit[2].symbol == :Y
        @test circuit[3].symbol == :Z
    end
    
    @testset "Two-Qubit Gates" begin
        qasm = """
        OPENQASM 2.0;
        qreg q[3];
        cx q[0],q[1];
        cnot q[1],q[2];
        cz q[0],q[2];
        """
        circuit = parse_qasm(qasm)
        
        @test length(circuit) == 3
        @test circuit[1] isa CliffordGate
        @test circuit[1].symbol == :CNOT
        @test circuit[1].qinds == [1, 2]
        @test circuit[2].symbol == :CNOT
        @test circuit[2].qinds == [2, 3]
        @test circuit[3].symbol == :CZ
        @test circuit[3].qinds == [1, 3]
    end
    
    @testset "SWAP Gate" begin
        qasm = """
        OPENQASM 2.0;
        qreg q[2];
        swap q[0],q[1];
        """
        circuit = parse_qasm(qasm)
        
        @test length(circuit) == 1
        @test circuit[1].symbol == :SWAP
        @test circuit[1].qinds == [1, 2]
    end
    
    @testset "Parametrized Rotation Gates" begin
        qasm = """
        OPENQASM 2.0;
        qreg q[3];
        rx(1.57) q[0];
        ry(3.14) q[1];
        rz(0.785) q[2];
        """
        circuit = parse_qasm(qasm)
        
        @test length(circuit) == 3
        @test circuit[1] isa FrozenGate
        @test circuit[1].gate isa PauliRotation
        @test circuit[1].gate.symbols == [:X]
        @test circuit[1].gate.qinds == [1]
        @test circuit[1].parameter ≈ 1.57
        
        @test circuit[2].gate.symbols == [:Y]
        @test circuit[2].gate.qinds == [2]
        @test circuit[2].parameter ≈ 3.14
        
        @test circuit[3].gate.symbols == [:Z]
        @test circuit[3].gate.qinds == [3]
        @test circuit[3].parameter ≈ 0.785
    end
    
    @testset "Parameter Expressions with π" begin
        qasm = """
        OPENQASM 2.0;
        qreg q[2];
        rx(pi/2) q[0];
        ry(pi) q[1];
        """
        circuit = parse_qasm(qasm)
        
        @test length(circuit) == 2
        @test circuit[1].parameter ≈ π/2
        @test circuit[2].parameter ≈ π
    end
    
    @testset "Skip Measurements and Barriers" begin
        qasm = """
        OPENQASM 2.0;
        qreg q[2];
        creg c[2];
        h q[0];
        barrier q;
        cx q[0],q[1];
        measure q[0] -> c[0];
        measure q[1] -> c[1];
        """
        circuit = parse_qasm(qasm)
        
        # Should only have h and cx gates, measurements ignored
        @test length(circuit) == 2
        @test circuit[1].symbol == :H
        @test circuit[2].symbol == :CNOT
    end
    
    @testset "Comments and Empty Lines" begin
        qasm = """
        OPENQASM 2.0;
        // This is a comment
        qreg q[2];
        
        h q[0];  // Hadamard on qubit 0
        
        cx q[0],q[1];
        // Another comment
        """
        circuit = parse_qasm(qasm)
        
        @test length(circuit) == 2
    end
    
    @testset "S Gate" begin
        qasm = """
        OPENQASM 2.0;
        qreg q[1];
        s q[0];
        """
        circuit = parse_qasm(qasm)
        
        @test length(circuit) == 1
        @test circuit[1].symbol == :S
    end
    
    @testset "Bell State Circuit" begin
        # Classic example: create Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2
        qasm = """
        OPENQASM 2.0;
        include "qelib1.inc";
        qreg q[2];
        creg c[2];
        h q[0];
        cx q[0],q[1];
        measure q[0] -> c[0];
        measure q[1] -> c[1];
        """
        circuit = parse_qasm(qasm)
        
        @test length(circuit) == 2
        @test circuit[1] isa CliffordGate
        @test circuit[1].symbol == :H
        @test circuit[2] isa CliffordGate
        @test circuit[2].symbol == :CNOT
    end
    
    @testset "GHZ State Circuit" begin
        # 3-qubit GHZ state
        qasm = """
        OPENQASM 2.0;
        qreg q[3];
        h q[0];
        cx q[0],q[1];
        cx q[1],q[2];
        """
        circuit = parse_qasm(qasm)
        
        @test length(circuit) == 3
        @test circuit[1].symbol == :H
        @test circuit[2].symbol == :CNOT
        @test circuit[3].symbol == :CNOT
    end
    
    @testset "parse_qasm_file" begin
        # Create a temporary QASM file
        tmpfile = tempname() * ".qasm"
        try
            qasm_content = """
            OPENQASM 2.0;
            qreg q[2];
            h q[0];
            cx q[0],q[1];
            """
            write(tmpfile, qasm_content)
            
            circuit = parse_qasm_file(tmpfile)
            
            @test length(circuit) == 2
            @test circuit[1].symbol == :H
            @test circuit[2].symbol == :CNOT
        finally
            # Clean up
            if isfile(tmpfile)
                rm(tmpfile)
            end
        end
    end
    
    @testset "Mixed Gate Types" begin
        qasm = """
        OPENQASM 2.0;
        qreg q[3];
        h q[0];
        rx(pi/4) q[0];
        cx q[0],q[1];
        ry(pi/2) q[1];
        cz q[1],q[2];
        z q[2];
        """
        circuit = parse_qasm(qasm)
        
        @test length(circuit) == 6
        @test circuit[1] isa CliffordGate
        @test circuit[2] isa FrozenGate
        @test circuit[3] isa CliffordGate
        @test circuit[4] isa FrozenGate
        @test circuit[5] isa CliffordGate
        @test circuit[6] isa CliffordGate
    end
    
end
