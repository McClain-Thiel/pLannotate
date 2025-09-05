# Publishing pLannotate to PyPI

This guide explains how to publish the pLannotate package to PyPI.

## Prerequisites

1. **PyPI Account**: Create accounts on both:
   - [PyPI](https://pypi.org/account/register/) (production)
   - [TestPyPI](https://test.pypi.org/account/register/) (testing)

2. **Install Publishing Tools**:
   ```bash
   pip install twine build
   # OR with uv (already available)
   ```

3. **API Tokens**: Set up API tokens for secure publishing:
   - Go to [PyPI Account Settings](https://pypi.org/manage/account/) → "API tokens"
   - Create a new token with scope "Entire account" 
   - Save the token securely

## Step-by-Step Publishing Process

### 1. Prepare the Package

Ensure all metadata is correct in `pyproject.toml`:
- ✅ Version number updated
- ✅ Description accurate
- ✅ GitHub URLs correct
- ✅ Dependencies specified with versions
- ✅ Classifiers appropriate

### 2. Test the Build

```bash
# Clean previous builds
rm -rf dist/ build/

# Build the package
uv build
# OR: python -m build

# Verify build contents
ls -la dist/
```

Should create:
- `plannotate-1.2.3.tar.gz` (source distribution)
- `plannotate-1.2.3-py3-none-any.whl` (wheel)

### 3. Test on TestPyPI First

```bash
# Upload to TestPyPI first
twine upload --repository testpypi dist/*

# Test installation from TestPyPI
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ plannotate

# Test that it works
python -c "from plannotate.annotate import annotate; print('✓ Package works')"
```

### 4. Publish to PyPI

If TestPyPI works correctly:

```bash
# Upload to real PyPI
twine upload dist/*

# Test installation from PyPI
pip install plannotate

# Verify it works
python -c "from plannotate.annotate import annotate; print('✓ Published successfully')"
```

## Authentication

### Option 1: API Tokens (Recommended)

Create a `.pypirc` file in your home directory:
```ini
[distutils]
index-servers = 
    pypi
    testpypi

[pypi]
repository = https://upload.pypi.org/legacy/
username = __token__
password = pypi-your-api-token-here

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = pypi-your-testpypi-token-here
```

### Option 2: Environment Variables

```bash
export TWINE_USERNAME=__token__
export TWINE_PASSWORD=pypi-your-api-token-here

# Then upload
twine upload dist/*
```

## Package Information

- **Package Name**: `plannotate`
- **Current Version**: `1.2.3`
- **Python Requirements**: `>=3.9`
- **Main Dependencies**: `biopython`, `pandas`, `reportlab`, `pyyaml`

## Post-Publishing Checklist

After successful PyPI publishing:

1. **Test Installation**:
   ```bash
   pip install plannotate
   python -c "from plannotate.annotate import annotate; print('Success!')"
   ```

2. **Update README**: Change installation instructions to use PyPI
3. **Create GitHub Release**: Tag the version and create a release
4. **Update Documentation**: Any references to installation methods
5. **Announce**: Share with the community

## Common Issues & Solutions

### Package Name Already Exists
- Check if `plannotate` is available on PyPI
- If taken, consider names like `plannotate-mcclain`, `py-plannotate`, etc.

### Large Package Size
- The package (~27MB) includes data files, which is normal for bioinformatics tools
- Consider if any large files can be moved to external downloads

### Missing Dependencies
- Make sure all imports work with specified dependency versions
- Test in a clean environment

### Metadata Issues
```bash
# Check package metadata
twine check dist/*
```

## Version Management

For future releases:

1. **Update Version**: Change in both `pyproject.toml` and `plannotate/__init__.py`
2. **Update Changelog**: Document changes
3. **Test Thoroughly**: Run full test suite
4. **Build & Upload**: Follow the same process

## Example Commands Summary

```bash
# Complete publishing workflow
rm -rf dist/
uv build
twine check dist/*
twine upload --repository testpypi dist/*
# Test installation & functionality
twine upload dist/*  # Upload to PyPI
```

## Resources

- [Python Packaging Guide](https://packaging.python.org/)
- [Twine Documentation](https://twine.readthedocs.io/)
- [PyPI Help](https://pypi.org/help/)
- [TestPyPI](https://test.pypi.org/)

---

**Note**: This package includes data files and external tool dependencies. Make sure users understand they need to install Diamond, BLAST, and Infernal separately for full functionality.