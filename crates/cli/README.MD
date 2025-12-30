# Install `b`

Pick a directory:

```bash
INSTALL_DIR="$HOME/.local/bin"
mkdir -p "$INSTALL_DIR"
```

---

## 1) Download

### macOS — arm64

```bash
curl -fL "https://raw.githubusercontent.com/josoriom/b/main/crates/cli/artifacts/macos-arm64" -o "$INSTALL_DIR/b"
chmod +x "$INSTALL_DIR/b"
```

### macOS — x86_64

```bash
curl -fL "https://raw.githubusercontent.com/josoriom/b/main/crates/cli/artifacts/macos-x86_64" -o "$INSTALL_DIR/b"
chmod +x "$INSTALL_DIR/b"
```

### Linux — x86_64

```bash
curl -fL "https://raw.githubusercontent.com/josoriom/b/main/crates/cli/artifacts/linux-x86_64" -o "$INSTALL_DIR/b"
chmod +x "$INSTALL_DIR/b"
```

### Linux — arm64

```bash
curl -fL "https://raw.githubusercontent.com/josoriom/b/main/crates/cli/artifacts/linux-arm64" -o "$INSTALL_DIR/b"
chmod +x "$INSTALL_DIR/b"
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
b --help
```

## 4) [DOCS](./src/README.MD)
