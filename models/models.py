from sqlalchemy import Column, Integer, String, LargeBinary, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

DATABASE_URL = 'postgresql://user:pass@db:5432/docking_db'
engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(bind=engine)

Base = declarative_base()

class DockingResult(Base):
    __tablename__ = 'docking_results'

    id = Column(Integer, primary_key=True, index=True)
    smiles = Column(String)
    protein_id = Column(String)
    result = Column(String)

def save_docking_result(smiles, protein_id, result):
    session = SessionLocal()
    docking_result = DockingResult(smiles=smiles, protein_id=protein_id, result=result)
    session.add(docking_result)
    session.commit()
    session.close()