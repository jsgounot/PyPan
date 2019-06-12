import argparse
import pypan

def launch(config_file, species, ncore) :
    config = config_file
    pypan.utils.set_config(config)
    
    fnames = pypan.nrfinder.nrclean(species, ncore=ncore)
    fnames = pypan.sample_sim.run(species, fnames, ncore=ncore)
    
    pypan.prediction.launch_fnames(species, fnames, ncore=ncore)
    pypan.soft_seg.run(species, fnames, ncore=ncore)
    pypan.sample_merge.run(species, fnames, ncore=ncore)

if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description='PyPan basic launcher')
    parser.add_argument("config", type=str, help="Path of your configuration file")
    parser.add_argument("species", type=str, help="Species to use")
    parser.add_argument("--ncore", type=int, help="Number of cores to use", default=1)

    args = parser.parse_args()

    launch(args.config, args.species, args.ncore)