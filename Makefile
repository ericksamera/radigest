BIN_DIR := bin
PREFIX ?= /usr/local
GO ?= go
PYTHON ?= python3

VERSION  := $(shell git describe --tags --dirty --always 2>/dev/null || echo dev)
GOFLAGS ?= -trimpath
LDFLAGS := -s -w -X main.version=$(VERSION)
SCRIPTS := scripts/radigest-screen-pairs scripts/radigest-rank-pairs scripts/radigest-plan-depth scripts/radigest-fit-size-model
CACHED_SCREEN_BIN := $(BIN_DIR)/radigest-screen-pairs-cached

.PHONY: all build install test lint tidy clean

all: build

build:
	mkdir -p $(BIN_DIR)
	$(GO) build $(GOFLAGS) -ldflags "$(LDFLAGS)" -o $(BIN_DIR)/radigest ./cmd/radigest
	$(GO) build $(GOFLAGS) -ldflags "$(LDFLAGS)" -o $(CACHED_SCREEN_BIN) ./cmd/radigest-screen-pairs-cached
	cp $(SCRIPTS) $(BIN_DIR)/
	chmod 0755 $(BIN_DIR)/radigest $(CACHED_SCREEN_BIN) $(BIN_DIR)/radigest-screen-pairs $(BIN_DIR)/radigest-rank-pairs $(BIN_DIR)/radigest-plan-depth $(BIN_DIR)/radigest-fit-size-model

install: build
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 0755 $(BIN_DIR)/radigest $(DESTDIR)$(PREFIX)/bin/radigest
	install -m 0755 $(CACHED_SCREEN_BIN) $(DESTDIR)$(PREFIX)/bin/radigest-screen-pairs-cached
	install -m 0755 $(BIN_DIR)/radigest-screen-pairs $(DESTDIR)$(PREFIX)/bin/radigest-screen-pairs
	install -m 0755 $(BIN_DIR)/radigest-rank-pairs $(DESTDIR)$(PREFIX)/bin/radigest-rank-pairs
	install -m 0755 $(BIN_DIR)/radigest-plan-depth $(DESTDIR)$(PREFIX)/bin/radigest-plan-depth
	install -m 0755 $(BIN_DIR)/radigest-fit-size-model $(DESTDIR)$(PREFIX)/bin/radigest-fit-size-model

test:
	$(GO) test $(GOFLAGS) ./... -count=1
	$(PYTHON) scripts/radigest_plan_depth_test.py

lint:
	@command -v golangci-lint >/dev/null 2>&1 || { \
	  echo "golangci-lint not found; install from https://golangci-lint.run/"; exit 0; }
	golangci-lint run ./...

tidy:
	$(GO) mod tidy

clean:
	rm -rf $(BIN_DIR)
