import configparser
import os.path

def read_config(dataset):
    assert os.path.exists('configs/' + dataset), "config file does not exist"
    config = configparser.ConfigParser()
    config.read('configs/' + dataset)

    config = config._sections
    config['PCA']['n_components'] = int(config['PCA']['n_components'])
    config['UMAP']['n_components'] = int(config['UMAP']['n_components'])
    config['UMAP']['n_neighbors'] = int(config['UMAP']['n_neighbors'])
    config['dataset']['w'] = int(config['dataset']['w'])
    config['dataset']['h'] = int(config['dataset']['h'])
    return config

def write_config(config_dict, filename):
    assert not os.path.exists('configs/' + filename), "config file exists"
    config = configparser.ConfigParser()
    config.read_dict(config_dict)
    with open('configs/' + filename, 'w') as fl:
        config.write(fl)