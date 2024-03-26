"""
    carbohydrate_parser.py
    Dylan Ross
    2019/07/17

        Utility for parsing carbohydrate names -- simple regex match for the patterns: 
            '[(]Hex[)][0-9]+'
            '^([LD][-])*[A-Za-z]+ose$'
            '^[A-Za-z]+itol$'
            '.*Cyclodextrin'
            'Lacto[-]N[-][A-Za-z]+ose.*'
"""


import re


def parse_carbohydrate(name):
    """
parse_carbohydrate
    description:
        Checks a compound name for matches with regex patterns for carbohydrates, returns
        a boolean reflecting whether the name looks like a carbohydrate.
    parameters:
        name (str) -- compound name
    returns:
        (bool) -- name matches carbohydrate regex pattern
"""
    if re.match(r'[(]Hex[)][0-9]+', name):
        return True

    if re.match(r'^([LD][-])*[A-Za-z]+ose$', name):
        return True

    if re.match(r'^[A-Za-z]+itol$', name):
        return True

    if re.match(r'.*Cyclodextrin', name):
        return True

    if re.match(r'Lacto[-]N[-][A-Za-z]+ose.*', name):
        return True
    
    return False
