-- C3SDB_schema.sql
-- Dylan H. Ross
-- 2019/01/07
--
--      defines the structures of the C3S.db Sqlite3 database


-- master table in which to combine all of the datasets
CREATE TABLE master (
    -- global unique string identifier
    g_id TEXT UNIQUE NOT NULL,
    -- compound name
    name TEXT NOT NULL,
    -- MS adduct
    adduct TEXT NOT NULL,
    -- mass, z, m/z and CCS
    mass REAL NOT NULL,
    z INTEGER NOT NULL,
    mz REAL NOT NULL,
    ccs REAL NOT NULL,
    -- neutral smiles structure
    smi TEXT,
    -- (rough) chemical class label
    chem_class_label TEXT,
    -- tag referencing which dataset the value is from
    src_tag TEXT NOT NULL,
    -- CCS type (DT, TW, ...)
    ccs_type TEXT NOT NULL,
    -- describe method used for CCS measurement (e.g. stepped-field, calibrated with polyalanine)
    ccs_method TEXT NOT NULL
);
