PLUGIN_NAME   = engine_physics
BUILD_DIR     = build
INSTALL_DIR   = $(HOME)/.hajimu/plugins/$(PLUGIN_NAME)
OUTPUT        = $(BUILD_DIR)/$(PLUGIN_NAME).hjp

NCPU = $(shell sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)

.PHONY: all clean install uninstall

all: $(OUTPUT)

$(OUTPUT): $(BUILD_DIR)/Makefile
	cmake --build $(BUILD_DIR) -j$(NCPU)

$(BUILD_DIR)/Makefile: CMakeLists.txt
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=Release -Wno-dev

install: all
	@mkdir -p $(INSTALL_DIR)
	cp $(OUTPUT) $(INSTALL_DIR)/
	@echo "  インストール: $(INSTALL_DIR)/$(PLUGIN_NAME).hjp"

uninstall:
	rm -f $(INSTALL_DIR)/$(PLUGIN_NAME).hjp

clean:
	rm -rf $(BUILD_DIR)
