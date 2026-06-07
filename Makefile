BIN_DIR := bin
PREFIX ?= /usr/local
GO ?= go
PYTHON ?= python3

VERSION  := $(shell git describe --tags --dirty --always 2>/dev/null || echo dev)
GOFLAGS ?= -trimpath
LDFLAGS := -s -w -X main.version=$(VERSION)
PUBLIC_SCRIPTS := scripts/radigest-fit-size-model
DEV_SCRIPTS := scripts/radigest-screen-pairs scripts/radigest-rank-pairs scripts/radigest-plan-depth
CACHED_SCREEN_BIN := $(BIN_DIR)/radigest-screen-pairs-cached
BENCH_SCREEN_BIN := $(BIN_DIR)/radigest-bench-screen-cached
DESIGN_BIN := $(BIN_DIR)/radigest-design

.PHONY: all build build-dev install install-dev test lint tidy clean

all: build

build:
	mkdir -p $(BIN_DIR)
	$(GO) build $(GOFLAGS) -ldflags "$(LDFLAGS)" -o $(BIN_DIR)/radigest ./cmd/radigest
	$(GO) build $(GOFLAGS) -ldflags "$(LDFLAGS)" -o $(DESIGN_BIN) ./cmd/radigest-design
	cp $(PUBLIC_SCRIPTS) $(BIN_DIR)/
	chmod 0755 $(BIN_DIR)/radigest $(DESIGN_BIN) $(BIN_DIR)/radigest-fit-size-model

build-dev: build
	$(GO) build $(GOFLAGS) -ldflags "$(LDFLAGS)" -o $(CACHED_SCREEN_BIN) ./cmd/radigest-screen-pairs-cached
	$(GO) build $(GOFLAGS) -ldflags "$(LDFLAGS)" -o $(BENCH_SCREEN_BIN) ./cmd/radigest-bench-screen-cached
	cp $(DEV_SCRIPTS) $(BIN_DIR)/
	chmod 0755 $(CACHED_SCREEN_BIN) $(BENCH_SCREEN_BIN) $(BIN_DIR)/radigest-screen-pairs $(BIN_DIR)/radigest-rank-pairs $(BIN_DIR)/radigest-plan-depth

install: build
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 0755 $(BIN_DIR)/radigest $(DESTDIR)$(PREFIX)/bin/radigest
	install -m 0755 $(DESIGN_BIN) $(DESTDIR)$(PREFIX)/bin/radigest-design
	install -m 0755 $(BIN_DIR)/radigest-fit-size-model $(DESTDIR)$(PREFIX)/bin/radigest-fit-size-model

install-dev: build-dev
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 0755 $(BIN_DIR)/radigest $(DESTDIR)$(PREFIX)/bin/radigest
	install -m 0755 $(DESIGN_BIN) $(DESTDIR)$(PREFIX)/bin/radigest-design
	install -m 0755 $(BIN_DIR)/radigest-fit-size-model $(DESTDIR)$(PREFIX)/bin/radigest-fit-size-model
	install -m 0755 $(CACHED_SCREEN_BIN) $(DESTDIR)$(PREFIX)/bin/radigest-screen-pairs-cached
	install -m 0755 $(BENCH_SCREEN_BIN) $(DESTDIR)$(PREFIX)/bin/radigest-bench-screen-cached
	install -m 0755 $(BIN_DIR)/radigest-screen-pairs $(DESTDIR)$(PREFIX)/bin/radigest-screen-pairs
	install -m 0755 $(BIN_DIR)/radigest-rank-pairs $(DESTDIR)$(PREFIX)/bin/radigest-rank-pairs
	install -m 0755 $(BIN_DIR)/radigest-plan-depth $(DESTDIR)$(PREFIX)/bin/radigest-plan-depth

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
