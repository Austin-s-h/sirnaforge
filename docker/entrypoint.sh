#!/bin/bash
# siRNAforge Docker entrypoint script
# Provides intelligent defaults and environment setup

set -euo pipefail

# Function to show available tools
show_environment() {
    echo "🧬 siRNAforge Comprehensive Environment"
    echo "======================================"
    echo ""
    echo "📦 Available tools:"
    echo "  • $(sirnaforge version 2>/dev/null || echo 'not available')"
    echo "  • nextflow: $(nextflow -version 2>/dev/null | head -n1 || echo 'not available')"
    echo "  • bwa-mem2: $(bwa-mem2 version 2>/dev/null | head -n1 || echo 'not available')"
    echo "  • samtools: $(samtools --version 2>/dev/null | head -n1 || echo 'not available')"
    echo "  • ViennaRNA: $(RNAfold --version 2>/dev/null | head -n1 || echo 'not available')"
    echo ""
    echo "🚀 Quick start examples:"
    echo "  # Design siRNAs for a gene"
    echo "  sirnaforge workflow TP53 --output-dir results"
    echo ""
    echo "  # Design from FASTA file"
    echo "  sirnaforge design input.fasta -o output.tsv"
    echo ""
    echo "  # Interactive shell"
    echo "  bash"
    echo ""
    echo "  # Get help for any command"
    echo "  sirnaforge --help"
    echo "  sirnaforge workflow --help"
}

# Function to check if we're in an interactive terminal
is_interactive() {
    [[ -t 0 && -t 1 && -t 2 ]]
}

# Main entrypoint logic
main() {
    # If no arguments provided, show smart default behavior
    if [[ $# -eq 0 ]]; then
        if is_interactive; then
            echo "🧬 Starting interactive siRNAforge environment..."
            show_environment
            echo ""
            echo "💡 You're in an interactive shell. Try 'sirnaforge --help' to get started."
            exec bash
        else
            # Non-interactive (e.g., in CI/CD), show environment info and exit
            show_environment
            return 0
        fi
    fi

    # If first argument is a known sirnaforge command, prepend 'sirnaforge'
    case "$1" in
        workflow|design|search|validate|config|version)
            exec sirnaforge "$@"
            ;;
        # If it's a shell command, execute directly
        bash|sh|/bin/bash|/bin/sh)
            exec "$@"
            ;;
        # If it's asking for help on the container
        help|--help|-h)
            show_environment
            echo ""
            echo "🔧 Container-specific commands:"
            echo "  help, --help     Show this help"
            echo "  bash, sh         Start interactive shell"
            echo "  env              Show environment info"
            echo ""
            sirnaforge --help
            ;;
        # Show environment info
        env|environment|info)
            show_environment
            ;;
        # Otherwise, execute the command as-is
        *)
            exec "$@"
            ;;
    esac
}

# Run the main function with all arguments
main "$@"
