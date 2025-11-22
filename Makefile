BIN_DIR := bin
GO ?= go

VERSION  := $(shell git describe --tags --dirty --always 2>/dev/null || echo dev)
GOFLAGS ?= -trimpath
LDFLAGS := -s -w -X main.version=$(VERSION)

.PHONY: all build test lint tidy clean

all: build

build:
	mkdir -p $(BIN_DIR)
	$(GO) build $(GOFLAGS) -ldflags "$(LDFLAGS)" -o $(BIN_DIR)/radigest ./cmd/radigest

test:
	$(GO) test $(GOFLAGS) ./... -count=1

lint:
	@command -v golangci-lint >/dev/null 2>&1 || { \
	  echo "golangci-lint not found; install from https://golangci-lint.run/"; exit 0; }
	golangci-lint run ./...

tidy:
	$(GO) mod tidy

clean:
	rm -rf $(BIN_DIR)
