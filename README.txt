Foreign Antigenic Space database query script & dependencies

Main script takes takes an input file (specified in script or as commandline argument) and queries the Schumi FAS database with the entries it parses

Required commandline arguments:
(1) input file path & filename (csv)
(2) tissue type
(3) lookup progress (specifies where to start in input file, starts at passed index +1; automatically obtained from FAS database through masterLookupController.R)
