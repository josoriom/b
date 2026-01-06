# Install `octo`

Pick a directory:

```bash
INSTALL_DIR="$HOME/.local/bin"
mkdir -p "$INSTALL_DIR"
```

---

## 1) Download

### macOS — arm64

```bash
curl -fL "https://raw.githubusercontent.com/josoriom/octo/main/crates/cli/artifacts/macos-arm64/octo" -o "$INSTALL_DIR/octo"
chmod +x "$INSTALL_DIR/octo"
```

### macOS — x86_64

```bash
curl -fL "https://raw.githubusercontent.com/josoriom/octo/main/crates/cli/artifacts/macos-x86_64/octo" -o "$INSTALL_DIR/octo"
chmod +x "$INSTALL_DIR/octo"
```

### Linux — x86_64

```bash
curl -fL "https://raw.githubusercontent.com/josoriom/octo/main/crates/cli/artifacts/linux-x86_64/octo" -o "$INSTALL_DIR/octo"
chmod +x "$INSTALL_DIR/octo"
```

### Linux — arm64

```bash
curl -fL "https://raw.githubusercontent.com/josoriom/octo/main/crates/cli/artifacts/linux-arm64/octo" -o "$INSTALL_DIR/octo"
chmod +x "$INSTALL_DIR/octo"
```

---

## 2) Add to PATH (.bashrc/.zshrc)

```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"
source "$HOME/.bashrc"
```

---

## 3) Test

```bash
octo --help
```

## 4) [DOCS](./src/README.MD)
