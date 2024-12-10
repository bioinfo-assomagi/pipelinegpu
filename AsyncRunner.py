import multiprocessing
class AsyncRunner():

    def __init__(self, store_results = True, **kwargs) -> None:
        self.samples = kwargs.pop("samples")
        self.store_results = store_results
        self.kwargs = kwargs
        # self.worker = worker
        
    def run(self, worker):
        q = self._generate_queue() if self.store_results else None
        args = self._prepare_args(q, self.kwargs)
        print(type(args))
        with multiprocessing.Pool(16) as p:
            p.map(worker, args)
            p.close()
            p.join()

        return q

    def _prepare_args(self, q, kwargs):
        # args = [(sample, q, kwargs) for sample in self.samples]
        args = []
        for sample in self.samples:
            kwargs.update({"sample": sample, "queue": q})
            args.append(kwargs)

        return args

    def set_store_results(self, store_results: bool):
        self.store_results = store_results

    def _generate_queue(self):
        return multiprocessing.Manager().Queue()