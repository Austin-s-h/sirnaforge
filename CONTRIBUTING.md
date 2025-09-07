# Contributing to siRNAforge

Thank you for your interest in contributing to siRNAforge! 🧬

## 🎯 Ways to Contribute

- 🐛 **Bug Reports**: Found a bug? Open an issue with details
- 💡 **Feature Requests**: Have an idea? We'd love to hear it
- 📖 **Documentation**: Help improve our docs
- 🧪 **Code**: Submit pull requests for bug fixes or features
- 🔬 **Testing**: Help test new features and report issues

## 🚀 Getting Started

### 1. Development Setup

```bash
# Fork the repository on GitHub, then clone
git clone https://github.com/YOUR-USERNAME/sirnaforge.git
cd sirnaforge

# Set up development environment
make install-dev
# or manually:
uv sync --dev
uv run pre-commit install
```

### 2. Development Workflow

```bash
# Create a feature branch
git checkout -b feature/your-feature-name

# Make your changes and add tests
# Run the development checks
make check  # Runs linting and tests

# Format your code
make format

# Commit with descriptive messages
git commit -m "feat: add new siRNA scoring algorithm"

# Push and create pull request
git push origin feature/your-feature-name
```

## 📋 Development Standards

### Code Quality
- **Type Hints**: All code must include proper type annotations
- **Testing**: New features require corresponding tests
- **Documentation**: Update docstrings and README as needed
- **Formatting**: Code is automatically formatted with Black and Ruff

### Testing
```bash
# Run all tests
make test

# Run tests with coverage
make test-cov

# Run specific test files
uv run pytest tests/unit/test_specific_module.py -v
```

### Linting & Formatting
```bash
# Check code quality
make lint

# Auto-format code
make format
```

## 📝 Pull Request Guidelines

### Before Submitting
- [ ] Code follows project style guidelines
- [ ] Tests pass locally (`make test`)
- [ ] New features include tests
- [ ] Documentation is updated
- [ ] Commit messages are descriptive

### PR Description Template
```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation update

## Testing
- [ ] Tests pass locally
- [ ] New tests added for new functionality

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
```

## 🐛 Bug Reports

When filing a bug report, please include:

1. **Environment details**:
   - Python version
   - Operating system
   - Package version

2. **Steps to reproduce**:
   - Minimal code example
   - Input files (if applicable)
   - Expected vs actual behavior

3. **Error messages**:
   - Full traceback
   - Log files if relevant

## 💡 Feature Requests

For new features, please:

1. **Check existing issues** to avoid duplicates
2. **Describe the use case** - why is this needed?
3. **Propose implementation** - how should it work?
4. **Consider scope** - should this be core functionality or a plugin?

## 🎨 Code Style

We use modern Python best practices:

- **Black** for code formatting (120 char line length)
- **Ruff** for linting and import sorting
- **MyPy** for type checking
- **Pre-commit hooks** for automated checks

## 🧪 Testing Guidelines

### Test Structure
- Unit tests in `tests/unit/`
- Integration tests in `tests/integration/`
- Pipeline tests in `tests/pipeline/`

### Test Coverage
- Aim for >90% code coverage
- Focus on critical paths and edge cases
- Include both positive and negative test cases

### Mock External Dependencies
- Mock file I/O operations
- Mock network requests
- Mock external tools (ViennaRNA, alignment tools)

## 📚 Documentation

- **Code**: Comprehensive docstrings with examples
- **API**: Auto-generated from docstrings
- **Tutorials**: Practical examples in `/examples`
- **README**: Keep the main README up-to-date

## 🚢 Release Process

1. Update `CHANGELOG.md`
2. Bump version in `src/sirnaforge/__init__.py`
3. Create GitHub release
4. CI automatically publishes to PyPI

## 🤝 Community Guidelines

- Be respectful and inclusive
- Help others learn and contribute
- Provide constructive feedback
- Follow the project's Code of Conduct

## 📞 Getting Help

- 💬 **GitHub Discussions**: For questions and ideas
- 🐛 **Issues**: For bugs and feature requests
- 📧 **Email**: For private matters

## 🏷️ Commit Message Format

We follow conventional commits:

- `feat:` - New features
- `fix:` - Bug fixes
- `docs:` - Documentation changes
- `style:` - Code style changes
- `refactor:` - Code refactoring
- `test:` - Adding or updating tests
- `chore:` - Maintenance tasks

Examples:
```
feat: add thermodynamic stability scoring
fix: handle edge case in GC content calculation
docs: update installation instructions
test: add integration tests for workflow
```

---

**Thank you for contributing to siRNAforge!** 🧬✨
