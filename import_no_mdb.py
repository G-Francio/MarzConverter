import classes as c
import defaults as d

user = c.QUser()


def get_observation_data(nameList):
    """
    Fallback method if mariadb is not found.

    Returns a list of neutral items in order to preserve the same script used with the DB itself.

    Parameters:
    - nameList (list): A list of names for which neutral observation data will be returned.

    Returns:
    list: A list containing neutral observation data for the given names, preserving compatibility with the database script.
    """
    return d.get_default_data(nameList)
