import warnings

import mariadb as mdb
import numpy as np

import classes as c
import defaults as d

user = c.QUser()


def authenticate_user(_user=user):
    """
    Authenticate a user for QubricsDB.

    Parameters:
    - _user: An object representing the user. Defaults to the global 'd.user'.

    Returns:
    tuple: A tuple containing the user object and the user's password.

    This function attempts to retrieve the username and password from the provided user object (_user).
    If the username or password is not available, it prompts the user to input the username and password
    using the 'input' and 'getpass' functions. The resulting tuple contains the user object and the password.
    """

    _uname, _pwd = _user.get_credentials()

    if _uname is None or _pwd is None:
        _user, _pwd = _user.get_credentials()
        # This will return None, None in case the user does not have access to the qubrics DB

    return _user, _pwd


# -------------------------------- ** -------------------------------- #


def DBConnect(user, db="QDB2"):
    """
    Connects to QubricsDB, returns a cursor and the connection object.
    """
    conn = mdb.connect(
        user=user.uname, password=user.pwd, host="127.0.0.1", port=3306, database=db
    )
    cur = conn.cursor()
    return cur, conn


# -------------------------------- ** -------------------------------- #


def get_observation_data(namelist, user=user, db="QDB2"):
    """
    Queries the DB for complementary information available per target.

    Can either find data based on the QID or, if observed by us, from the named spectrum.

    Parameters:
    - namelist (list): A list of targets. Targets can be specified either by QID or by named spectrum.

    Returns:
    list: A list containing complementary observation data for the specified targets. The data is queried from the database.
    If no data is found based on QID, it attempts to retrieve data based on the named spectrum.

    Note:
    - QID (Query ID): A numerical identifier used for specifying targets in the database.
    - The function internally distinguishes targets based on whether they are specified by QID or named spectrum.
    """

    _name_list_QIDs = [
        e for e in namelist if (e.isdecimal() or e.split("_")[0].isdecimal())
    ]
    _name_list_files = [e for e in namelist if not e.isdecimal()]

    data = _get_data_from_QID(_name_list_QIDs, user, db=db)
    if len(data) == 0:
        return _get_observation_data(_name_list_files, user, db=db)

    return data


# -------------------------------- ** -------------------------------- #


def _get_observation_data(name_list, user=user, db="QDB2"):
    """
    Gathers data for observations that needs conversion.
    `nameList` is a list of fileNames taken from Drive/QSOCandidates.
    """
    observation_data = []
    file_name_list = []

    try:
        cur, conn = DBConnect(user, db=db)
        cur.execute("SELECT file FROM Qubrics.Observations")
        for c in cur:
            file_name_list.append(c)
        file_name_list = np.array(file_name_list)[:, 0]

        for name in name_list:
            idxs = np.where(np.char.find(file_name_list, "/" + name) > 0)

            if len(idxs[0]) > 1:
                print("Multiple matches for {}, using fallback data.".format(name))
                observation_data.append(d.get_default_data_single(name))
            elif len(idxs[0]) == 0:
                print(
                    "Can not find spec data on DB for {}, fallback on mock data.".format(
                        name
                    )
                )
                observation_data.append(d.get_default_data_single(name))
            else:
                fileName = file_name_list[idxs[0][0]]
                cur.execute(
                    "SELECT target_qid, RAd, DECd, otypeid, z_spec, targetflag, qflag, notes FROM Qubrics.Observations WHERE file = '{}'".format(
                        fileName
                    )
                )
                for c in cur:
                    c_objtype = (
                        c[0],  # qid
                        c[1],  # ra
                        c[2],  # dec
                        d.typedict(c[3]),  # otypeid -> string id
                        c[4],  # z_spec
                        c[5],  # targetflag
                        c[6],  # qflag
                        c[7],  # notes
                    )
                    observation_data.append(c_objtype)

        conn.close()
        return np.array(observation_data)

    except mdb.OperationalError:
        print("Connection to the DB failed, mock data will be used")
        return d.get_default_data(name_list)


# -------------------------------- ** -------------------------------- #


def _get_data_from_QID(qid_list, user=user, db="QDB2"):
    """
    Queries the DB if the name is numerical (and thus a QID).

    This allows downloading complementary information even if the target is not from
    our observations but from literature.

    Parameters:
    - qid_list (list): A list of numerical names (QIDs) to query from the database.
    - user: An object representing the user for database connection.
    - db (str): The database name. Defaults to "QDB2".

    Returns:
    numpy.ndarray: An array containing complementary observation data for the specified QIDs.
    The data includes QID, RAd, DECd, otypeid, and z_spec.

    Note:
    - QID (Query ID): A numerical identifier used for specifying targets in the database.
    - The function retrieves information from the 'Qubrics.All_info' table in the specified database.
    - In case of unexpected query results or errors, the function issues warnings and uses fallback data.
    """
    if len(qid_list) == 0:
        return np.array([])

    try:
        qid_list_out = []
        cur, conn = DBConnect(user, db=db)

        for _qid in qid_list:
            # Can't query all info at once, in the unlikely case that a numerical name is not in the DB.
            # Also, this will return wrong info in that case. Can't do anything about it, need to check
            #  beforehand!

            cur.execute(
                "SELECT qid, RAd, DECd, otypeid, z_spec FROM Qubrics.All_info WHERE qid = ?",
                (_qid,),
            )
            query_res = [*cur]
            # Check if I get exactly one result out of the query, otherwise warn the user and fallback on
            #  mock data
            if len(query_res) == 1:
                qid_list_out.append(np.array([*query_res[0]] + ["", "A", ""]))
            elif len(query_res) == 0:
                warnings.warn(
                    f"Query did not return additional information for qid: {_qid}"
                )
                qid_list_out.append(d.get_default_data_single(_qid))
            elif len(query_res) > 1:
                warnings.warn(
                    f"Too many results for qid: {_qid}, this should never happen!\nUsing fallback data."
                )
                qid_list_out.append(d.get_default_data_single(_qid))
            conn.close()

        return np.array(qid_list_out)

    except mdb.OperationalError:
        print("Connection to the DB failed, mock data will be used")
        return d.get_default_data(qid_list)
