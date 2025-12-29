# Publishing ROSE2 to PyPI

## Pre-Publication Checklist

- [x] Package structure created with `pyproject.toml`
- [x] All code modernized with type hints
- [x] Tests created and passing (20/20 tests pass)
- [x] README updated with installation instructions
- [x] Version set to 2.0.0
- [x] Package builds successfully
- [x] CLI commands work correctly

## Build the Package

The package has already been built. The distribution files are in `dist/`:

```bash
ls -lh dist/
# rose2-2.0.0-py3-none-any.whl (16M)
# rose2-2.0.0.tar.gz (16M)
```

To rebuild:

```bash
# Clean previous builds
rm -rf dist/ build/ *.egg-info

# Build new distributions
python3 -m build
```

## Test the Package Locally

### Install in Development Mode

```bash
pip install -e .
```

### Test Installation from Wheel

```bash
pip install dist/rose2-2.0.0-py3-none-any.whl
```

### Run Tests

```bash
pytest tests/ -v
```

All 20 tests should pass.

### Test CLI Commands

```bash
# Test main CLI
rose2 --version
rose2 --help

# Test subcommands
rose2 main --help
rose2 bamToGFF --help
rose2 geneMapper --help

# Test standalone entry points
rose2-main --help
rose2-bamToGFF --help
rose2-geneMapper --help
```

## Publish to PyPI

### 1. Create PyPI Account

- Sign up at https://pypi.org/account/register/
- Enable 2FA (required)
- Generate an API token at https://pypi.org/manage/account/token/

### 2. Create ~/.pypirc

```ini
[pypi]
username = __token__
password = pypi-your-api-token-here
```

### 3. Upload to Test PyPI (Recommended First)

```bash
# Upload to Test PyPI
python3 -m twine upload --repository testpypi dist/*

# Test install from Test PyPI
pip install --index-url https://test.pypi.org/simple/ rose2
```

### 4. Upload to Production PyPI

```bash
# Check the package first
python3 -m twine check dist/*

# Upload to PyPI
python3 -m twine upload dist/*
```

### 5. Verify Publication

```bash
# Install from PyPI
pip install rose2

# Test it works
rose2 --version
```

## Post-Publication

### Create GitHub Release

1. Tag the release:
```bash
git tag -a v2.0.0 -m "Release version 2.0.0 - Modernized Python package"
git push origin v2.0.0
```

2. Create release on GitHub:
- Go to https://github.com/stjude/ROSE2/releases/new
- Select tag `v2.0.0`
- Title: "ROSE2 v2.0.0 - Modernized PyPI Package"
- Description: Include changelog and installation instructions

### Update Documentation

- [ ] Update repository README with PyPI badge
- [ ] Add link to PyPI package page
- [ ] Update installation instructions
- [ ] Add citation information

### Announce

Consider announcing on:
- Bioinformatics mailing lists
- Twitter/X
- Biostars
- Reddit (r/bioinformatics)

## Version Updates

For future releases:

1. Update version in `pyproject.toml`
2. Update `__version__` in `rose2/__init__.py`
3. Update CHANGELOG/README
4. Build and test
5. Publish to PyPI
6. Create GitHub release

## Maintenance

### Monitoring

- Watch PyPI download statistics: https://pypistats.org/packages/rose2
- Monitor GitHub issues: https://github.com/stjude/ROSE2/issues
- Check for security vulnerabilities

### Updates

Keep dependencies updated:
```bash
pip install --upgrade pip setuptools wheel build twine
```

## Troubleshooting

### Common Issues

**Problem**: Upload fails with "File already exists"
**Solution**: Increment version number in `pyproject.toml` and rebuild

**Problem**: Import errors after installation
**Solution**: Check that annotation files are included in MANIFEST.in

**Problem**: CLI commands not found
**Solution**: Verify entry points are correct in `pyproject.toml`

## Contact

For publishing questions:
- Ming Tang: tangming2005@gmail.com
- GitHub Issues: https://github.com/stjude/ROSE2/issues
