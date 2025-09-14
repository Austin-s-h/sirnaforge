#!/bin/bash
# siRNAforge Local Development Setup Script
# Sets up a complete conda environment for local development

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

# Check if conda/mamba/micromamba is available
check_conda() {
    if command -v micromamba &> /dev/null; then
        CONDA_CMD="micromamba"
        print_info "Using micromamba (fastest option)"
    elif command -v mamba &> /dev/null; then
        CONDA_CMD="mamba"
        print_info "Using mamba (faster than conda)"
    elif command -v conda &> /dev/null; then
        CONDA_CMD="conda"
        print_info "Using conda"
    else
        print_error "Neither conda, mamba, nor micromamba found. Please install one of:"
        echo "  • micromamba (recommended): https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html"
        echo "  • Mambaforge: https://mamba.readthedocs.io/en/latest/installation.html"
        echo "  • Miniconda: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
}

# Main setup function
setup_environment() {
    print_info "Setting up siRNAforge local development environment..."

    # Check conda availability
    check_conda

    # Create environment
    print_info "Creating conda environment 'sirnaforge-dev'..."
    $CONDA_CMD env create -f environment-dev.yml

    # Activate environment and install the package
    print_info "Activating environment and installing siRNAforge..."
    eval "$($CONDA_CMD shell.bash hook)"
    $CONDA_CMD activate sirnaforge-dev

    # Install the package itself (since conda environment has dependencies but not the package)
    print_info "Installing siRNAforge package..."
    pip install -e .

    print_success "Environment setup complete!"
    print_info "To activate the environment in future sessions:"
    echo "  conda activate sirnaforge-dev"
    echo ""
    print_info "To verify installation:"
    echo "  sirnaforge version"
    echo "  sirnaforge config"
    echo ""
    print_info "To run tests:"
    echo "  make test-local-python  # Fast Python-only tests"
    echo "  make test               # Full test suite (may take longer)"
}

# Function to update environment
update_environment() {
    print_info "Updating siRNAforge development environment..."

    check_conda

    eval "$($CONDA_CMD shell.bash hook)"
    $CONDA_CMD activate sirnaforge-dev

    print_info "Updating conda environment..."
    $CONDA_CMD env update -f environment-dev.yml

    print_info "Reinstalling siRNAforge package..."
    pip install -e .

    print_success "Environment updated!"
}

# Function to clean up
cleanup_environment() {
    print_warning "This will remove the sirnaforge-dev conda environment."
    read -p "Are you sure? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        check_conda
        $CONDA_CMD env remove -n sirnaforge-dev
        print_success "Environment removed."
    else
        print_info "Operation cancelled."
    fi
}

# Help function
show_help() {
    echo "siRNAforge Local Development Environment Setup"
    echo ""
    echo "Usage: $0 [COMMAND]"
    echo ""
    echo "Commands:"
    echo "  setup     Create and configure the development environment (default)"
    echo "  update    Update existing environment with latest dependencies"
    echo "  cleanup   Remove the development environment"
    echo "  help      Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0                    # Setup environment"
    echo "  $0 setup             # Same as above"
    echo "  $0 update            # Update environment"
    echo "  $0 cleanup           # Remove environment"
}

# Main script logic
case "${1:-setup}" in
    "setup")
        setup_environment
        ;;
    "update")
        update_environment
        ;;
    "cleanup")
        cleanup_environment
        ;;
    "help"|"-h"|"--help")
        show_help
        ;;
    *)
        print_error "Unknown command: $1"
        echo ""
        show_help
        exit 1
        ;;
esac
