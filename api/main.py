from fastapi import FastAPI, BackgroundTasks, HTTPException
from pydantic import BaseModel
from celery.result import AsyncResult
from workers.tasks import dock_compound_task, celery_app
from models.molecule import MoleculeProcessor

# Configure Celery with Redis backend
celery_app.conf.update(
    result_backend='redis://redis:6379/0'
)

app = FastAPI()

class DockingJobRequest(BaseModel):
    smiles: str
    protein_id: str

@app.post("/submit_job")
async def submit_job(job_request: DockingJobRequest):
    # Validate SMILES and generate 3D structure before submitting to worker
    mol = MoleculeProcessor.smiles_to_3d(job_request.smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string or 3D structure generation failed")
        
    task = dock_compound_task.delay(job_request.smiles, job_request.protein_id)
    return {"task_id": task.id}

@app.get("/job_status/{task_id}")
async def job_status(task_id: str):
    result = AsyncResult(task_id)
    if result.state == 'PENDING':
        return {"status": "Pending"}
    elif result.state == 'SUCCESS':
        return {"status": "Completed", "result": result.result}
    elif result.state == 'FAILURE':
        return {"status": "Failed", "error": str(result.info)}
    else:
        return {"status": result.state}