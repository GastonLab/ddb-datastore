import ConfigParser
from collections import defaultdict


def configure_runtime(infile):
    """Parse the configuration settings from a file
    :param infile: input filename
    :type infile: string.
    :returns:  dict -- A configuration dictionary.
    """

    configuration = defaultdict()
    config = ConfigParser.SafeConfigParser()
    config.read(infile)

    try:
        config.options('settings')
    except ConfigParser.NoSectionError:
        sys.stderr.write("No section: settings in file\n")
        sys.exit()

    try:
        config.options('resources')
    except ConfigParser.NoSectionError:
        sys.stderr.write("No section: resources in file\n")
        sys.exit()

    # Set all options specified in file
    for option in config.options('settings'):
            configuration[option] = config.get('settings', option)

    for resource in config.options('resources'):
            configuration[resource] = config.get('resources', resource)

    # Configure individual tools
    for section in config.sections():
        if section != 'settings' and section != 'resources':
            tool = section
            options = config.options(tool)
            tool_dict = dict()

            # Set all specified options
            for option in options:
                tool_dict[option] = config.get(tool, option)

            configuration[tool] = tool_dict

    return configuration


def configure_samples(infile, configuration):
    """Parse the sample-level configuration settings from a file
    :param infile: input filename
    :type infile: string.
    :param configuration: project/run level configuration
    :type configuration: dictionary
    :returns:  dict -- A configuration dictionary.
    """

    samples = dict()

    config = ConfigParser.SafeConfigParser()
    config.read(infile)

    for sample in config.sections():
        sample_dict = dict()
        for option in config.options(sample):
            sample_dict[option] = config.get(sample, option)

        if 'regions' not in sample_dict.keys():
            sample_dict['regions'] = configuration['regions']
        if 'snv_regions' not in sample_dict.keys():
            if 'snv_regions' in configuration:
                sample_dict['snv_regions'] = configuration['snv_regions']
        if 'indel_regions' not in sample_dict.keys():
            if 'indel_regions' in configuration:
                sample_dict['indel_regions'] = configuration['indel_regions']

        samples[sample] = sample_dict

    return samples