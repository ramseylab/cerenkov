from urllib.parse import urlencode
from urllib.request import Request, urlopen
import pandas
from io import StringIO


class SnapQuery:
    # 1000 Genomes Pilot 1 / HapMap release 22 / HapMap release 21
    _dataset_options = {'onekgpilot', 'rel22', 'rel21'}

    _population_options = {'onekgpilot': {'CEU', 'YRI', 'CHBJPT'},
                           'rel22': {'CEU', 'YRI', 'JPT+CHB'},
                           'rel21': {'CEU', 'YRI', 'JPT+CHB'}}

    @property
    def dataset_options(self):
        return type(self)._dataset_options

    @property
    def population_options(self):
        return type(self)._population_options

    def __init__(self):
        self.snp_dataset = None
        self.population_panels = None
        self.r_squared = None
        self.distance_limit = None
        self.snp_list = None

    def dataset(self, _dataset):
        if _dataset not in self.dataset_options:
            raise ValueError("[SnapQuery] Please choose dataset among {}".format(self.dataset_options))
        self.snp_dataset = _dataset
        return self

    def population(self, _population):
        _population_set = set(_population)
        if any([_p not in self.population_options[self.snp_dataset] for _p in _population_set]):
            raise ValueError("[SnapQuery] Please choose population among {} when using {} dataset".
                             format(self.population_options[self.snp_dataset], self.snp_dataset))
        self.population_panels = _population_set
        return self

    def r_squared_threshold(self, _value):
        self.r_squared = _value
        return self

    def distance_limit_in_kb(self, _value):
        self.distance_limit = _value
        return self

    def snp(self, _snp_seq):
        self.snp_list = ",".join(_snp_seq)
        return self

    def _execute(self, population, verbose=False):
        url = 'http://archive.broadinstitute.org/mpg/snap/ldsearch.php'
        param = {
            'snpList': self.snp_list,
            'hapMapRelease': self.snp_dataset,
            'hapMapPanel': population,
            'RSquaredLimit': self.r_squared,
            'distanceLimit': self.distance_limit * 1000,
            'downloadType': 'file',
            # 'includeQuerySnp': 'on',
            'arrayFilter': 'query',
            'submit': 'search',
            'suppressWarnings': 'on',
            'columnList[]': 'GP',
        }

        # if verbose:
        #     print("[SnapQuery] dataset = '{dataset}', population = '{population}', r_squared = {r2}, "
        #           "distance = {distance}kb, snp_list = '{snp}'".format(dataset=self.snp_dataset,
        #                                                                population=population,
        #                                                                r2=self.r_squared,
        #                                                                distance=self.distance_limit,
        #                                                                snp=self.snp_list))

        data = urlencode(param)
        data = data.encode('ascii')  # data should be bytes
        req = Request(url, data)

        with urlopen(req) as response:
            page = response.read()
            dfm = pandas.read_csv(StringIO(page.decode("utf8")), encoding='utf8', header=0, sep="\t",
                                  skipinitialspace=True)

            if verbose:
                print("[SnapQuery] dataset = '{dataset}', population = '{population}', "
                      "r_squared = {r2}, distance = {distance}kb, #input_entries = {n_in}, "
                      "#output_entries = {n_out}".format(dataset=self.snp_dataset,
                                                         population=population,
                                                         r2=self.r_squared,
                                                         distance=self.distance_limit,
                                                         n_in=len(self.snp_list),
                                                         n_out=dfm.shape[0]))
            return dfm

    def execute(self, verbose=False):
        dfm_list = [self._execute(population, verbose) for population in self.population_panels]
        if len(dfm_list) == 1:
            return dfm_list[0]
        else:
            dfm = pandas.concat(dfm_list, ignore_index=True, axis=0)
            return dfm
