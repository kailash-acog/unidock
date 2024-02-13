IMAGE_NAME?="dockerhub.aganitha.ai:4443/dev/unidock-image-kailash"
CONTAINER_NAME?=$(USER)-unidock

build-image:
	@cd devops && docker build -t $(IMAGE_NAME) .

start-container:
	@cd devops && IMAGE_NAME=$(IMAGE_NAME) CONTAINER_NAME=$(CONTAINER_NAME) \
		./start_container.sh;

enter-container:
	@echo "You are inside the Container: \033[1;33m$(CONTAINER_NAME)\033[0m"; \
	docker exec -u root -it $(CONTAINER_NAME) bash || true


## PostgreSQL 
PG_IMAGE_NAME?=unidock-postgres-kailash
PG_CONTAINER_NAME?=${USER}-unidock-postgres

enter-container-pg:
	@echo "You are inside the Container: \033[1;33m$(PG_CONTAINER_NAME)\033[0m"; \
	docker exec -u root -it $(PG_CONTAINER_NAME) bash || true

enter-db:
	@docker exec -it $(PG_CONTAINER_NAME) psql -U demo_user -d demo_db || true

