#!/bin/bash

set -e

INSTALL_DIR="$HOME/.local/bin"
REPO_URL="https://github.com/RG8420/probfun.git"
TEMP_DIR=$(mktemp -d)

echo "=== Probability Calculator (probcalc) Installer ==="

if [[ "$1" == "--uninstall" ]]; then
    echo "[*] Uninstalling probcalc..."
    rm -f "$INSTALL_DIR/probcalc"
    echo "  Uninstalled."
    exit 0
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
if [[ -f "$SCRIPT_DIR/src/main.c" && -f "$SCRIPT_DIR/Makefile" ]]; then
    echo "[*] Running from local build..."
    BUILD_DIR="$SCRIPT_DIR"
else
    BUILD_DIR="$TEMP_DIR/probcalc"
    echo "[*] Cloning repository..."
    git clone --depth 1 "$REPO_URL" "$BUILD_DIR"
fi

cd "$BUILD_DIR"

echo "[1/3] Building probcalc..."
make clean
make

echo "[2/3] Installing probcalc..."
mkdir -p "$INSTALL_DIR"
cp probcalc "$INSTALL_DIR/"

echo "[3/3] Configuring PATH..."

if [[ ":$PATH:" != *":$INSTALL_DIR:"* ]]; then
    if [[ -f "$HOME/.bashrc" ]]; then
        if ! grep -q "# probcalc" "$HOME/.bashrc" 2>/dev/null; then
            echo "" >> "$HOME/.bashrc"
            echo "# probcalc" >> "$HOME/.bashrc"
            echo "export PATH=\"\$HOME/.local/bin:\$PATH\"" >> "$HOME/.bashrc"
        fi
    fi
    if [[ -f "$HOME/.profile" ]]; then
        if ! grep -q "# probcalc" "$HOME/.profile" 2>/dev/null; then
            echo "" >> "$HOME/.profile"
            echo "# probcalc" >> "$HOME/.profile"
            echo "export PATH=\"\$HOME/.local/bin:\$PATH\"" >> "$HOME/.profile"
        fi
    fi
    echo "  Added \$HOME/.local/bin to PATH"
fi

echo ""
echo "  Installation successful!"
echo ""
echo "  Run 'probcalc' to start!"
echo ""
echo "  To uninstall, run:  bash install.sh --uninstall"

rm -rf "$TEMP_DIR"