function [y,x] = get_utah_array_coordinates_from_channel_ID(channel_ID)

    electrode_map_matrix =  [0 45 88 84 80 76 72 68 60 0;
                            47 6  2  8  86 82 78 64 56 52;
                            41 43 21 23 4  66 74 70 71 67;
                            39 37 17 19 7  11 62 50 69 63;
                            35 33 13 15 91 3  54 58 65 59;
                            31 29 9  5  87 83 95 57 61 55;
                            25 27 1  93 85 79 75 96 92 53;
                            46 48 44 89 81 77 73 94 90 49;
                            42 40 36 32 28 24 20 16 12 51;
                            0  38 34 30 26 22 18 14 10 0 ];
    
   [y,x] = find(electrode_map_matrix == channel_ID);
end