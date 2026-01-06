CARGO_BIN         := octo
CLI_MANIFEST      := Cargo.toml
ARTIFACTS         := artifacts

WORKSPACE_ROOT    := $(abspath $(CURDIR)/../..)
WORKSPACE_TARGET  := $(WORKSPACE_ROOT)/target

DOCKER_IMAGE      := rust:1-bullseye

WINDOWS_ARM64_IMAGE  := dockcross/windows-arm64:latest
WINDOWS_ARM64_TARGET := aarch64-pc-windows-gnullvm

.PHONY: all mac linux windows macos-arm64 macos-x86_64 linux-amd64 linux-arm64 windows-amd64 windows-arm64 clean

all: mac linux windows

mac: macos-arm64 macos-x86_64
linux: linux-amd64 linux-arm64
windows: windows-amd64 windows-arm64

macos-arm64:
	rustup target add aarch64-apple-darwin
	cargo build --manifest-path $(CLI_MANIFEST) --release --target aarch64-apple-darwin
	mkdir -p $(ARTIFACTS)/macos-arm64
	cp $(WORKSPACE_TARGET)/aarch64-apple-darwin/release/$(CARGO_BIN) $(ARTIFACTS)/macos-arm64/$(CARGO_BIN)

macos-x86_64:
	rustup target add x86_64-apple-darwin
	cargo build --manifest-path $(CLI_MANIFEST) --release --target x86_64-apple-darwin
	mkdir -p $(ARTIFACTS)/macos-x86_64
	cp $(WORKSPACE_TARGET)/x86_64-apple-darwin/release/$(CARGO_BIN) $(ARTIFACTS)/macos-x86_64/$(CARGO_BIN)

linux-amd64:
	docker run --rm --platform=linux/amd64 \
	  -e CARGO_TARGET_DIR=/work/target-linux-amd64 \
	  -v $$PWD/../..:/work -w /work/crates/cli \
	  --entrypoint /usr/local/cargo/bin/cargo $(DOCKER_IMAGE) \
	  build --manifest-path Cargo.toml --release
	mkdir -p $(ARTIFACTS)/linux-x86_64
	cp $(WORKSPACE_ROOT)/target-linux-amd64/release/$(CARGO_BIN) $(ARTIFACTS)/linux-x86_64/$(CARGO_BIN)

linux-arm64:
	docker run --rm --platform=linux/arm64 \
	  -e CARGO_TARGET_DIR=/work/target-linux-arm64 \
	  -v $$PWD/../..:/work -w /work/crates/cli \
	  --entrypoint /usr/local/cargo/bin/cargo $(DOCKER_IMAGE) \
	  build --manifest-path Cargo.toml --release
	mkdir -p $(ARTIFACTS)/linux-arm64
	cp $(WORKSPACE_ROOT)/target-linux-arm64/release/$(CARGO_BIN) $(ARTIFACTS)/linux-arm64/$(CARGO_BIN)

windows-amd64:
	docker run --rm --platform=linux/amd64 \
	  -e CARGO_TARGET_DIR=/work/target-windows-amd64 \
	  -v $$PWD/../..:/work -w /work/crates/cli $(DOCKER_IMAGE) bash -lc '\
	    set -euo pipefail; \
	    export PATH=/usr/local/cargo/bin:$$PATH; \
	    apt-get update && apt-get install -y --no-install-recommends \
	      curl ca-certificates gcc-mingw-w64-x86-64 g++-mingw-w64-x86-64 mingw-w64; \
	    command -v rustup >/dev/null 2>&1 || { curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal; . "$$HOME/.cargo/env"; }; \
	    rustup target add x86_64-pc-windows-gnu; \
	    RUSTFLAGS="-C linker=x86_64-w64-mingw32-gcc" \
	      cargo build --manifest-path Cargo.toml --release --target x86_64-pc-windows-gnu'
	mkdir -p $(ARTIFACTS)/windows-x86_64
	cp $(WORKSPACE_ROOT)/target-windows-amd64/x86_64-pc-windows-gnu/release/$(CARGO_BIN).exe $(ARTIFACTS)/windows-x86_64/$(CARGO_BIN).exe

windows-arm64:
	docker run --rm --platform=linux/amd64 \
	  -e CARGO_TARGET_DIR=/work/target-windows-arm64 \
	  -v $$PWD/../..:/work -w /work/crates/cli $(WINDOWS_ARM64_IMAGE) \
	  bash -lc 'set -euo pipefail; \
	    export PATH="$$HOME/.cargo/bin:/usr/local/cargo/bin:$$PATH"; \
	    command -v rustup >/dev/null 2>&1 || { curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal; . "$$HOME/.cargo/env"; }; \
	    rustup target add $(WINDOWS_ARM64_TARGET); \
	    cargo build --manifest-path Cargo.toml --release --target $(WINDOWS_ARM64_TARGET)'
	mkdir -p $(ARTIFACTS)/windows-arm64
	cp $(WORKSPACE_ROOT)/target-windows-arm64/$(WINDOWS_ARM64_TARGET)/release/$(CARGO_BIN).exe $(ARTIFACTS)/windows-arm64/$(CARGO_BIN).exe

clean:
	cargo clean --manifest-path $(CLI_MANIFEST)
	rm -rf $(ARTIFACTS) \
	  $(WORKSPACE_ROOT)/target-linux-amd64 \
	  $(WORKSPACE_ROOT)/target-linux-arm64 \
	  $(WORKSPACE_ROOT)/target-windows-amd64 \
	  $(WORKSPACE_ROOT)/target-windows-arm64
