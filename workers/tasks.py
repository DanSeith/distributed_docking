from celery import Celery
from docking.docking import perform_docking
from cache.cache import get_prepared_protein
from models.models import save_docking_result

celery_app = Celery('tasks', broker='redis://redis:6379/0')

@celery_app.task(bind=True, max_retries=3)
def dock_compound_task(self, smiles, protein_id):
    try:
        prepared_protein = get_prepared_protein(protein_id)
        docking_result = perform_docking(smiles, prepared_protein)
        save_docking_result(smiles, protein_id, docking_result)
        return docking_result
    except Exception as exc:
        raise self.retry(exc=exc, countdown=60)