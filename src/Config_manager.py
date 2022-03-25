import json


class Config_manager:
    def _get_data(self, key: str):
        with open("./src/config.json", "r") as f:
            return json.load(f)[key]

    def get_outputPath(self):
        return self._get_data("outputPath")
