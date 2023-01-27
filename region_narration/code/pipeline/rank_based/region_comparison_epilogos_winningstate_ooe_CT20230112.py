import numpy as np
import sys
import h5py


def main(mapped_path, expected_path, num_states, out_path):
    # load mapped winning state data
    f = open(mapped_path, 'r')
    winning_state_data = []
    for line in f.readlines():
        winning_state_data.append(
            np.float_(line.strip().split('|')[1].split(';'))
        )
    winning_state_data = np.array(winning_state_data)
    f.close()

    # count each state winning in each region
    counts = []
    for region in range(len(winning_state_data)):
        region_states = (winning_state_data[region][
            ~np.isnan(winning_state_data[region])
        ]).astype(int)
        state_count = np.bincount(region_states)
        state_count = np.append(state_count,
                                np.int_([0] * (num_states - len(state_count))))
        counts.append(state_count)
    counts = np.array(counts)

    # observed ratio
    observed = (counts.T / np.sum(counts, axis=1)).T

    # expected ratio (whole genome)
    f = open(expected_path, 'r')
    expected = np.float_(f.readlines()[2].strip().split('\t'))
    f.close()

    ooe = observed / expected

    # save ooe array
    ooe_f = h5py.File(out_path, 'w')
    ooe_f.create_dataset('base_arr', data=ooe)
    ooe_f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
