#!/usr/bin/env python

"""Rolling in the deep.

Example:

    ./batch_run.py -n 16 path/to/jobs.txt

The `-n` parameter specifies the number workers to use. This should be roughly similar to the number
of cores in your machine.

The `jobs.txt` file should contain a list of obj to process, one per line. The obj files should be
relative to the folder containing `jobs.txt`.
"""

import argparse
import asyncio
from pathlib import Path
import shutil

project_root = Path.cwd().parent

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("job_file", help="Job file")
    parser.add_argument(
        "-n", "--num-workers", type=int, default=4, help="Number of workers"
    )
    return parser.parse_args()


async def run_subprocess(command):
    """Run a subprocess command asynchronously and return the result."""
    proc = await asyncio.create_subprocess_shell(
        command,
        stderr=asyncio.subprocess.PIPE,
        stdout=asyncio.subprocess.PIPE,
        cwd=project_root / "build",
    )

    stdout, stderr = await proc.communicate()

    if proc.returncode == 0:
        pass
    else:
        print(f"[ERROR] {command}\n{stderr.decode()}")


async def worker(name, queue):
    """Worker that takes jobs from the queue and runs them."""
    while True:
        # Get a job (command) from the queue
        file = await queue.get()
        if file is None:
            # Exit signal received
            break
        print(f"Worker {name} is executing: {file}")
        executable = "./scriptIPC"
        num_samples = 1000
        ipc_repo = project_root.parent / "rigid-ipc"
        assert ipc_repo.exists() and ipc_repo.is_dir()
        assert (ipc_repo / "build" / "rigid_ipc_sim").exists()
        command = f"{executable} --samples {num_samples} --ipc --ipc_dir {ipc_repo} -m {file}"
        await run_subprocess(command)
        # Signal that the job is done
        queue.task_done()

        print(f"Worker {name} finished: {file}")


async def process_files(files: list[Path], num_workers: int = 4):
    print(f"Processing {len(files)} files with {num_workers} workers")

    # Create a queue with the commands to execute
    queue: asyncio.Queue = asyncio.Queue()

    # Enqueue all the commands into the queue
    for file in files:
        await queue.put(file)

    # Create and start worker tasks
    workers = [
        asyncio.create_task(worker(f"Worker-{i}", queue)) for i in range(num_workers)
    ]

    # Wait for the queue to be processed
    await queue.join()

    # Send exit signals to workers
    for _ in workers:
        await queue.put(None)

    # Wait for all workers to finish
    await asyncio.gather(*workers)


def main():
    args = parse_args()
    job_file = Path(args.job_file)
    assert job_file.exists()
    job_folder = job_file.parent
    files = []
    with open(job_file, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            obj_file = job_folder / line.strip()
            assert obj_file.exists()
            files.append(obj_file)
            logs_folder = obj_file.parent / (obj_file.stem + "_logs")
            if logs_folder.exists() and logs_folder.is_dir():
                backup_folder = logs_folder.parent / (logs_folder.name + "_old")
                if backup_folder.exists():
                    shutil.rmtree(backup_folder)
                shutil.move(logs_folder, backup_folder)
    asyncio.run(process_files(files, args.num_workers))


if __name__ == "__main__":
    main()