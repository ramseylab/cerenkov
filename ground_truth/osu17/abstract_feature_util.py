class AbstractFeatureUtil(object):
    def __init__(self):
        self._temp_dest = None
        self._db_config_key = None
        self._src_data_dir = None
        self._src_data_fn = None

    @property
    def temp_dest(self):
        return self._temp_dest

    @temp_dest.setter
    def temp_dest(self, value):
        self._temp_dest = value

    @property
    def db_config_key(self):
        return self._db_config_key

    @db_config_key.setter
    def db_config_key(self, value):
        self._db_config_key = value

    @property
    def src_data_dir(self):
        return self._src_data_dir

    @src_data_dir.setter
    def src_data_dir(self, value):
        self._src_data_dir = value

    @property
    def src_data_fn(self):
        return self._src_data_fn

    @src_data_fn.setter
    def src_data_fn(self, value):
        self._src_data_fn = value

    def get_feat(self, _input):
        raise NotImplementedError("This is an abstract method.")

    def save_temp(self, _result):
        raise NotImplementedError("This is an abstract method.")

    def extract(self, _input):
        _result = self.get_feat(_input)

        if self.temp_dest is not None:
            self.save_temp(_result)

        return _result
