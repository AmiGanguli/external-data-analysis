import random
import asyncio
import sys
import aiohttp
import logging


class Connection:
    def __init__(self, session, service_endpoint):
        self.service_endpoint = service_endpoint
        self.failure_count = 0
        self.session = session

    # TODO: Need to simulate failures somehow to test this.
    async def get(self, url, max_retry=10, min_wait=3):
        success = False
        retries = 0
        wait_time = min_wait
        response_json = None
        while not success:
            try:
                logging.debug(f'Starting GET request for {url}')
                async with self.session.get(url) as response:
                    response_json = await response.json()
                success = True
            except aiohttp.ClientConnectorError:
                # Connection errors are typically local network errors.
                # Just wait for one minute and then retry. Print a
                # notification rather than giving up.
                #
                logging.warn(
                    f"Connection error for {url}. Will re-try forever. Sleeping for 1 minute.",
                )
                await asyncio.sleep(60)
            except aiohttp.ClientConnectionError as e:
                # All other errors, we keep retrying for a long time.
                #
                if retries > max_retry:
                    logging.error(
                        f"Maximum number of attempts exceeded for {url}. Exit.",
                    )
                    raise SystemExit(e)
                logging.warn(
                    f"Warn: Request error for {url}. Retry #{retries} of {max_retry}. Sleeping for {wait_time} seconds.",
                )
                await asyncio.sleep(wait_time)
                retries += 1
                wait_time += random.randint(1, wait_time)
        logging.debug(f'Completed GET request for {url}')
        return response_json

    async def post(self, url, data, max_retry=10, min_wait=3):
        success = False
        retries = 0
        wait_time = min_wait
        response_json = None
        while not success:
            try:
                logging.debug(f'Starting POST request for {url}')
                async with self.session.post(url, data=data) as response:
                    response_json = await response.json()
                success = True
            except aiohttp.ClientConnectorError:
                # Connection errors are typically local network errors.
                # Just wait for one minute and then retry. Print a
                # notification rather than giving up.
                #
                print(
                    f"Warn: Connection error for {url}. Will re-try forever. Sleeping for 1 minute.",
                    file=sys.stderr,
                )
                await asyncio.sleep(60)
            except aiohttp.ClientConnectionError as e:
                # All other errors, we keep retrying for a long time.
                #
                if retries > max_retry:
                    print(
                        f"Error: Maximum number of attempts exceeded for {url}. Exit.",
                        file=sys.stderr,
                    )
                    raise SystemExit(e)
                print(
                    f"Warn: Request error for {url}. Retry #{retries} of {max_retry}. Sleeping for {wait_time} seconds.",
                    file=sys.stderr,
                )
                await asyncio.sleep(wait_time)
                retries += 1
                wait_time += random.randint(1, wait_time)
        return response_json

    async def getEventsHierarchy(self, tax_id):
        return await self.get(
            f'{self.service_endpoint}/data/eventsHierarchy/{tax_id}')

    async def getProductData(self, product_id):
        return await self.get(
            f'{self.service_endpoint}/data/query/{product_id}')

    async def getProductDataMultiple(self, product_ids):
        product_data = []
        product_ids = [str(id) for id in product_ids]
        for i in range(0, len(product_ids), 20):
            request_data = ','.join(product_ids[i:i+20])
            logging.debug(
                f'Request  {len(product_ids)} products: {request_data}.')
            result = await self.post(
                f'{self.service_endpoint}/data/query/ids',
                request_data,
            )
            logging.debug(f'product data post returned {result}')
            for item in result:
                product_data.append(item)
        return product_data

    async def getSpecies(self):
        return await self.get(f'{self.service_endpoint}/data/species/all')

    async def getParticipantsReferenceEntities(self, id):
        return await self.get(f'{self.service_endpoint}/data/participants/{id}/referenceEntities')

    async def getParticipantsPhysicalEntities(self, id):
        return await self.get(f'{self.service_endpoint}/data/participants/{id}/participatingPhysicalEntities')
