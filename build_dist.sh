python -m conda_index repo
sed "3i\ \ - file://$(pwd)/repo" environment_source.yml > environment.yml
voici build --XeusAddon.empack_config=empack_config.yaml --contents content --output-dir dist
# jupyter lite build --contents content --XeusAddon.empack_config=empack_config.yaml --output-dir dist